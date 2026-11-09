#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script de procesamiento quimioinformático completo:
Curado, cálculo de descriptores, huellas dactilares Morgan (ECFP4 - FP_2048),
análisis estadístico (Mann-Whitney U + Delta de Cliff) y visualización del
espacio químico mediante PCA y UMAP.

Versión optimizada y corregida.
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, QED, rdMolDescriptors, rdFingerprintGenerator

# Silenciar advertencias excesivas de RDKit para mantener limpia la terminal
RDLogger.DisableLog('rdApp.*')

try:
    import umap.umap_ as umap
except ImportError:
    try:
        import umap
    except ImportError:
        umap = None

# ==========================================
# 1. DEFINICIÓN DE FUNCIONES Y DICCIONARIOS
# ==========================================

DESCRIPTORES = {
    "MW":             Descriptors.MolWt,
    "cLogP":          Descriptors.MolLogP,
    "TPSA":           Descriptors.TPSA,
    "HBD":            Descriptors.NumHDonors,
    "HBA":            Descriptors.NumHAcceptors,
    "RotB":           Descriptors.NumRotatableBonds,
    "FraccionCsp3":   Descriptors.FractionCSP3,
    "Anillos":        rdMolDescriptors.CalcNumRings,
    "Estereocentros": lambda m: len(Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False)),
    "AtomosPesados":  Descriptors.HeavyAtomCount,
    "QED":            QED.qed,
}

# Generador ECFP4 (Morgan Radius 2, 2048 bits)
generador_morgan = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

def calcular_descriptores(smiles):
    """Calcula el conjunto de descriptores fisicoquímicos para una estructura SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {nombre: funcion(mol) for nombre, funcion in DESCRIPTORES.items()}

def huella_morgan(smiles):
    """Genera el vector Numpy de la huella dactilar Morgan ECFP4 (FP 2048, radio 2)."""
    mol = Chem.MolFromSmiles(smiles)
    return np.array(generador_morgan.GetFingerprint(mol)) if mol else None

def delta_cliff(a, b):
    """Proporción de pares en que a supera a b, menos la inversa."""
    a, b = np.asarray(a), np.asarray(b)
    mayor = sum((x > b).sum() for x in a)
    menor = sum((x < b).sum() for x in a)
    return (mayor - menor) / (len(a) * len(b))

def magnitud_cliff(d):
    """Determina la magnitud del tamaño de efecto Delta de Cliff."""
    d = abs(d)
    if d < 0.147: return "despreciable"
    if d < 0.33:  return "pequeña"
    if d < 0.474: return "mediana"
    return "grande"

# ==========================================
# 2. EJECUCIÓN PRINCIPAL DEL FLUJO DE TRABAJO
# ==========================================

def main():
    # --- Rutas basadas en el directorio actual del script o estructura limpia ---
    # Si prefieres rutas absolutas específicas, puedes reasignarlas aquí.
    base_dir = os.path.dirname(os.path.abspath(__file__)) if '__file__' in locals() else "."
    
    # Rutas relativas ajustables a tu estructura de directorios
    ruta_corpus = os.path.join(base_dir, "data/processed/torresprocessed/torresF2/curado.csv")
    ruta_aprobados = os.path.join(base_dir, "data/processed/torresprocessed/torresF2/referencia_aprobados.csv")
    ruta_aleatorios = os.path.join(base_dir, "data/processed/torresprocessed/torresF2/referencia_aleatoria.csv")
    
    dir_salida_datos = os.path.join(base_dir, "data/processed/torresprocessed/torresF2")
    dir_salida_resultados = os.path.join(base_dir, "results/torresresults/torresF2")
    dir_figures = os.path.join(dir_salida_resultados, "figures")
    
    ruta_salida_descriptores = os.path.join(dir_salida_datos, "descriptores.csv")
    ruta_res_principal = os.path.join(dir_salida_resultados, "comparacion_descriptores.csv")
    ruta_res_excluidos = os.path.join(dir_salida_resultados, "descriptores_excluidos.csv")
    
    ruta_fig_pca = os.path.join(dir_figures, "espacio_quimico_pca.png")
    ruta_fig_umap = os.path.join(dir_figures, "espacio_quimico_umap.png")

    os.makedirs(dir_salida_datos, exist_ok=True)
    os.makedirs(dir_salida_resultados, exist_ok=True)
    os.makedirs(dir_figures, exist_ok=True)

    # --- 1. Carga y consolidación de datos ---
    print("==================================================")
    print("1. CARGA Y CONSOLIDACIÓN DE DATOS")
    print("==================================================")
    corpus = pd.read_csv(ruta_corpus)
    aprobados = pd.read_csv(ruta_aprobados)
    aleatorios = pd.read_csv(ruta_aleatorios)

    corpus["conjunto"] = "productos_naturales"
    if "conjunto" not in aprobados.columns:
        aprobados["conjunto"] = "aprobados"
    if "conjunto" not in aleatorios.columns:
        aleatorios["conjunto"] = "aleatorios"

    todos = pd.concat([corpus, aprobados, aleatorios], ignore_index=True)

    cols_smiles = [c for c in todos.columns if 'smiles' in c.lower()]
    todos['SMILES_UNIFICADO'] = None
    for col in cols_smiles:
        todos['SMILES_UNIFICADO'] = todos['SMILES_UNIFICADO'].fillna(todos[col])

    print("Conteo de moléculas por conjunto:")
    print(todos["conjunto"].value_counts())
    print(f"Moléculas con SMILES válido: {todos['SMILES_UNIFICADO'].notna().sum()} / {len(todos)}\n")

    # --- 2. Cálculo de descriptores fisicoquímicos ---
    print("==================================================")
    print("2. CÁLCULO DE DESCRIPTORES FISICOQUÍMICOS")
    print("==================================================")
    calculados = todos['SMILES_UNIFICADO'].apply(lambda x: calcular_descriptores(x) if isinstance(x, str) else None)
    lista_dict = [res if isinstance(res, dict) else {} for res in calculados]
    df_descriptores = pd.DataFrame(lista_dict)

    for col in df_descriptores.columns:
        todos[col] = df_descriptores[col]

    print("--- Conteo de datos válidos por grupo y descriptor ---")
    print(todos.groupby("conjunto")[list(DESCRIPTORES.keys())].count())
    print()

    # --- 3. Cálculo de Huellas Dactilares ECFP4 (Morgan FP_2048) ---
    print("==================================================")
    print("3. CÁLCULO DE HUELLAS DACTILARES ECFP4 (MORGAN FP_2048)")
    print("==================================================")
    col_smiles = 'SMILES_UNIFICADO' if 'SMILES_UNIFICADO' in todos.columns else ('SMILES' if 'SMILES' in todos.columns else 'smiles')
    todos['Morgan_FP_2048'] = todos[col_smiles].apply(lambda x: huella_morgan(x) if isinstance(x, str) else None)

    print("Muestra de Huellas Morgan (ECFP4) para Productos Naturales (primeras 5 filas):")
    print(todos[todos['conjunto'] == 'productos_naturales'][[col_smiles, 'Morgan_FP_2048']].head())
    print()

    # --- 4. Exportar entregable de corpus curado con descriptores (Optimizado sin recálculo) ---
    print("==================================================")
    print("4. EXPORTACIÓN DE CORPUS CURADO Y DESCRIPTORES")
    print("==================================================")
    entregable_final = todos[todos['conjunto'] == 'productos_naturales'].copy()
    entregable_final.to_csv(ruta_salida_descriptores, index=False)
    print(f"✓ Archivo CSV exportado en disco:\n  {ruta_salida_descriptores}")
    print(f"  - Registros: {entregable_final.shape[0]} | Columnas: {entregable_final.shape[1]}\n")

    # --- 5. Análisis Estadístico ---
    print("==================================================")
    print("5. ANÁLISIS ESTADÍSTICO (MANN-WHITNEY U & DELTA DE CLIFF)")
    print("==================================================")
    filas = []
    excluidos = []

    for descriptor in DESCRIPTORES:
        pn = todos.loc[todos.conjunto == "productos_naturales", descriptor].dropna()
        if len(pn) < 2:
            excluidos.append({
                "descriptor": descriptor,
                "grupo": "productos_naturales",
                "n_observaciones": len(pn),
                "motivo": "Muestra insuficiente"
            })
            continue

        mediana_pn = pn.median()
        iqr_pn = pn.quantile(.75) - pn.quantile(.25)

        for referencia in ("aprobados", "aleatorios"):
            ref = todos.loc[todos.conjunto == referencia, descriptor].dropna()
            if len(ref) < 2:
                excluidos.append({
                    "descriptor": descriptor,
                    "grupo": referencia,
                    "n_observaciones": len(ref),
                    "motivo": "Muestra insuficiente"
                })
                continue

            u_stat, p = mannwhitneyu(pn, ref, alternative="two-sided")
            d = (2 * u_stat) / (len(pn) * len(ref)) - 1

            filas.append({
                "descriptor": descriptor,
                "contra": referencia,
                "mediana_pn": mediana_pn,
                "IQR_pn": iqr_pn,
                "mediana_ref": ref.median(),
                "IQR_ref": ref.quantile(.75) - ref.quantile(.25),
                "delta_cliff": round(d, 3),
                "magnitud": magnitud_cliff(d),
                "p": p,
            })

    tabla = pd.DataFrame(filas)
    if not tabla.empty:
        tabla["p_corregido"] = multipletests(tabla["p"], method="holm")[1]
        tabla.sort_values("delta_cliff", key=abs, ascending=False).to_csv(ruta_res_principal, index=False)
        print(f"✓ Comparación estadística guardada en disco ({len(tabla)} pares):\n  {ruta_res_principal}")

    tabla_excluidos = pd.DataFrame(excluidos).drop_duplicates(subset=["descriptor", "grupo"])
    tabla_excluidos.to_csv(ruta_res_excluidos, index=False)
    print(f"✓ Registro de excluidos guardado en disco:\n  {ruta_res_excluidos}\n")

    # --- 6. PCA (Figura 1: espacio_quimico_pca.png) ---
    print("==================================================")
    print("6. ANÁLISIS DE COMPONENTES PRINCIPALES (PCA)")
    print("==================================================")
    df_valid = todos.dropna(subset=list(DESCRIPTORES.keys())).copy()
    X_pca = df_valid[list(DESCRIPTORES.keys())]
    X_pca_escalado = StandardScaler().fit_transform(X_pca)

    pca = PCA(n_components=2, random_state=42)
    coords_pca = pca.fit_transform(X_pca_escalado)
    var = pca.explained_variance_ratio_ * 100
    print(f"Varianza explicada: PC1 {var[0]:.1f}% | PC2 {var[1]:.1f}%")

    cargas = pd.DataFrame(pca.components_.T, index=list(DESCRIPTORES.keys()), columns=["PC1", "PC2"])
    print("\n--- Cargas de los Componentes (Loadings) ---")
    print(cargas.round(2).sort_values("PC1", key=abs, ascending=False))

    # Generación de espacio_quimico_pca.png
    df_plot_pca = pd.DataFrame(coords_pca, columns=["PC1", "PC2"])
    df_plot_pca["conjunto"] = df_valid["conjunto"].values

    plt.figure(figsize=(9, 6))
    sns.scatterplot(
        data=df_plot_pca,
        x="PC1",
        y="PC2",
        hue="conjunto",
        style="conjunto",
        palette="Set2",
        alpha=0.7,
        s=60
    )
    plt.title("Espacio Químico (PCA): Productos Naturales vs Referencias", fontsize=12, fontweight='bold')
    plt.xlabel(f"PC1 ({var[0]:.1f}% varianza explicada)")
    plt.ylabel(f"PC2 ({var[1]:.1f}% varianza explicada)")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.legend(title="Conjunto de datos")
    plt.tight_layout()

    plt.savefig(ruta_fig_pca, dpi=300)
    plt.close()
    print(f"\n✓ Figura PCA guardada en disco:\n  {ruta_fig_pca}\n")

    # --- 7. UMAP sobre Huellas ECFP4 (Figura 2: espacio_quimico_umap.png) ---
    print("==================================================")
    print("7. VISUALIZACIÓN UMAP SOBRE HUELLAS DACTILARES ECFP4")
    print("==================================================")
    # Comprobación segura de arrays numpy para evitar ValueError en Pandas
    mask_fp = todos['Morgan_FP_2048'].apply(lambda x: isinstance(x, np.ndarray))
    df_fp_valid = todos[mask_fp].copy()
    
    if len(df_fp_valid) > 0:
        X_fp = np.stack(df_fp_valid['Morgan_FP_2048'].values)
        
        if umap is not None:
            reducer = umap.UMAP(n_components=2, metric='jaccard', random_state=42, n_neighbors=15, min_dist=0.1)
            coords_umap = reducer.fit_transform(X_fp)
        else:
            print("(!) Módulo umap-learn no instalado localmente en esta ejecución. Se utiliza TSNE como respaldo.")
            from sklearn.manifold import TSNE
            reducer = TSNE(n_components=2, random_state=42, metric='jaccard')
            coords_umap = reducer.fit_transform(X_fp)

        df_plot_umap = pd.DataFrame(coords_umap, columns=["UMAP 1", "UMAP 2"])
        df_plot_umap["conjunto"] = df_fp_valid["conjunto"].values

        plt.figure(figsize=(9, 6))
        sns.scatterplot(
            data=df_plot_umap,
            x="UMAP 1",
            y="UMAP 2",
            hue="conjunto",
            style="conjunto",
            palette="Set2",
            alpha=0.7,
            s=60
        )
        plt.title("Espacio Químico (UMAP sobre huellas ECFP4)", fontsize=12, fontweight='bold')
        plt.xlabel("UMAP 1")
        plt.ylabel("UMAP 2")
        plt.grid(True, linestyle="--", alpha=0.5)
        plt.legend(title="Conjunto de datos")
        plt.tight_layout()

        plt.savefig(ruta_fig_umap, dpi=300)
        plt.close()
        print(f"✓ Figura UMAP guardada en disco:\n  {ruta_fig_umap}")
    else:
        print("(!) No se encontraron huellas dactilares válidas para UMAP.")

    print("==================================================")
    print("PROCESAMIENTO FINALIZADO CON ÉXITO")
    print("==================================================")

if __name__ == "__main__":
    main()