#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import sys
import argparse
from collections import defaultdict

def parse_sam_to_matrix(sam_file, output_tsv):
    arrays = defaultdict(dict)
    max_spacer_idx = 0

    print(f"--- Procesando archivo SAM: {sam_file} ---")

    with open(sam_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue

            # Extraemos las variables sin usar corchetes (a prueba de terminal)
            qname = parts.pop(0)
            flag_str = parts.pop(0)
            ref_genome = parts.pop(0)

            flag = int(flag_str)

            # Ignoramos líneas de alineamientos secundarios/suplementarios 
            if flag & 256 or flag & 2048:
                continue

            if ref_genome == '*':
                final_target = "Unknown"
            else:
                targets = set()
                targets.add(ref_genome)

                # 1. Buscamos el NM (número de mutaciones) del alineamiento principal
                primary_nm = None
                for tag in parts:
                    if tag.startswith("NM:i:"):
                        nm_str = tag.replace("NM:i:", "")
                        primary_nm = int(nm_str)
                        break

                # 2. Buscamos alineamientos alternativos y evaluamos si son "igual de buenos"
                for tag in parts:
                    if tag.startswith("XA:Z:"):
                        # Formato XA: chr,pos,CIGAR,NM;
                        xa_content = tag.replace("XA:Z:", "")
                        alt_alignments = xa_content.split(';')
                        
                        for alt in alt_alignments:
                            if alt:
                                alt_parts = alt.split(',')
                                # Comprobamos que tenga las 4 partes esperadas
                                if len(alt_parts) >= 4:
                                    # Sacamos los datos secuencialmente sin usar corchetes
                                    alt_ref = alt_parts.pop(0)
                                    alt_pos = alt_parts.pop(0)
                                    alt_cigar = alt_parts.pop(0)
                                    alt_nm_str = alt_parts.pop(0)
                                    
                                    try:
                                        alt_nm = int(alt_nm_str)
                                        # CONDICIÓN CLAVE: Solo lo guardamos si tiene los mismos (o menos) errores
                                        if primary_nm is None or alt_nm <= primary_nm:
                                            targets.add(alt_ref)
                                    except ValueError:
                                        pass
                
                # Unimos todas las dianas válidas separadas por /
                final_target = " / ".join(sorted(targets))

            try:
                array_id, spacer_pos_str = qname.rsplit('_', 1)
                spacer_pos = int(spacer_pos_str)
            except ValueError:
                continue

            # Guardamos la diana
            current_array_dict = arrays.setdefault(array_id, dict())
            current_array_dict.update({spacer_pos: final_target})

            if spacer_pos > max_spacer_idx:
                max_spacer_idx = spacer_pos

    if not arrays:
        print("No se encontraron alineamientos válidos para procesar.")
        return

    print(f"  Generando matriz de salida: {output_tsv}")
    with open(output_tsv, 'w', newline='') as out_f:
        writer = csv.writer(out_f, delimiter='\t')

        # Cabecera con la columna Array_Length
        header = list()
        header.append('Array_ID')
        header.append('Array_Length')
        for i in range(1, max_spacer_idx + 1):
            header.append(f'Spacer_{i}')
        writer.writerow(header)

        def sort_key(k):
            try:
                return int(k)
            except ValueError:
                return k

        for array_id in sorted(arrays.keys(), key=sort_key):
            row = list()
            row.append(array_id)

            # Calculamos la longitud real
            array_dict = arrays.get(array_id)
            array_length = max(array_dict.keys())
            row.append(str(array_length))

            # Rellenamos las columnas
            for i in range(1, max_spacer_idx + 1):
                if i <= array_length:
                    row.append(array_dict.get(i, "Unknown"))
                else:
                    row.append("-")

            writer.writerow(row)

    print("--- Conversión completada con éxito. ---")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convierte SAM a TSV cruzando solo dianas de igual calidad.")
    parser.add_argument("sam_file", help="Archivo SAM de entrada.")
    parser.add_argument("-o", "--output", default="crispr_targets_matrix.tsv", help="Archivo TSV de salida.")

    args = parser.parse_args()
    parse_sam_to_matrix(args.sam_file, args.output)
