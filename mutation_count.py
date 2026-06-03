import collections
from Levenshtein import editops, distance

# --- CONFIGURACIÓN ---
spacer_ref = "CCGATGGGATCCGACGTTCGCCGTCTGCGTTACC"
archivo_fastq = "8S3ZPF_1_769_PCR2_2.fastq"
max_errores = 3
# ---------------------

conteo_variantes = collections.Counter()
stats_mutaciones = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0}
total_spacers = 0

print(f"Analizando: {spacer_ref}\n")

with open(archivo_fastq, "r") as f:
    for i, line in enumerate(f):
        if i % 4 == 1:
            read = line.strip().upper()
            l_ref = len(spacer_ref)
            
            mejor_dist = max_errores + 1
            mejor_subseq = ""
            
            # Buscamos la mejor coincidencia en la lectura
            for j in range(len(read) - l_ref + 1):
                subseq = read[j:j+l_ref]
                d = distance(spacer_ref, subseq)
                if d < mejor_dist:
                    mejor_dist = d
                    mejor_subseq = subseq
            
            if mejor_dist <= max_errores:
                conteo_variantes[mejor_subseq] += 1
                stats_mutaciones[mejor_dist] += 1
                total_spacers += 1

# 1. Tabla de variantes detalladas
print(f"{'Secuencia':<35} | {'Detalle Mutación':<30} | {'Cant.'}")
print("-" * 80)
for seq, count in conteo_variantes.most_common():
    ops = editops(spacer_ref, seq)
    detalles = [f"{op[0].upper()}:{spacer_ref[op[1]] if op[0]!='insert' else ''}->{seq[op[2]] if op[0]!='delete' else ''}(p{op[1]+1})" for op in ops]
    desc = ", ".join(detalles) if detalles else "PERFECTA"
    print(f"{seq:<35} | {desc:<30} | {count}")

# 2. RESUMEN POR NÚMERO DE MUTACIONES
print("\n" + "="*40)
print("   RESUMEN POR NÚMERO DE MUTACIONES")
print("="*40)
for n_mut in range(max_errores + 1):
    cantidad = stats_mutaciones[n_mut]
    porcentaje = (cantidad / total_spacers * 100) if total_spacers > 0 else 0
    etiqueta = "Perfectos" if n_mut == 0 else f"{n_mut} Mutación/es"
    print(f"{etiqueta:<15}: {cantidad:>6} ({porcentaje:>5.1f}%)")

print("-" * 40)
print(f"SUMA TOTAL ENCONTRADA: {total_spacers}")
print("="*40)
