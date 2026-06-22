# Guide de démarrage rapide

## 1. Préparer le dossier du projet

Décompressez le dossier `thalamic-diaschisis-pipeline-v1.0.0` dans un emplacement distinct des données patients.

Ne placez jamais les DICOM, NIfTI, tableaux cliniques ou résultats patients dans le dépôt GitHub.

## 2. Créer la configuration locale

Copiez :

```text
config/example_config.m
```

et renommez la copie :

```text
config/local_config.m
```

Dans `local_config.m`, adaptez les quatre chemins obligatoires :

```matlab
cfg.root_dir      = 'D:\Vincent\CT_perf';
cfg.spm_dir       = 'C:\...\spm12';
cfg.ctseg_dir     = 'C:\...\CTseg';
cfg.fieldtrip_dir = 'C:\...\fieldtrip';
```

Vérifiez que `cfg.ctseg_dir` contient ou permet de trouver :

```text
spm_CTseg.m
mu_CTseg.nii
```

## 3. Vérifier les données d’entrée

Chaque patient doit avoir une structure similaire à :

```text
D:\Vincent\CT_perf\
├── Sub001\
│   ├── ANAT\
│   └── Perf_T\
├── Sub002\
│   ├── ANAT\
│   └── Perf_T\
└── ...
```

Les dossiers `ANAT` et `Perf_T` contiennent les fichiers DICOM correspondants.

## 4. Lancer le pipeline complet

Dans MATLAB :

```matlab
cd('CHEMIN\VERS\thalamic-diaschisis-pipeline-v1.0.0')
run_all
```

Le pipeline réalise successivement :

1. la conversion DICOM vers NIfTI ;
2. le recalage et la normalisation ;
3. la création des masques et cartes de groupe ;
4. l’extraction des valeurs Tmax thalamiques ;
5. des contrôles automatiques de cohérence.

## 5. Lancer uniquement une partie

Conversion et normalisation :

```matlab
PerfCT_Reg_MNI_001
```

Extraction thalamique :

```matlab
Extract_ROI_Thalamus_001
```

## 6. Vérifier les résultats

Ouvrez :

```text
D:\Vincent\CT_perf\Group_results\
```

Vérifiez notamment :

```text
dicom_conversion_summary.csv
processing_summary.csv
thalamus_CT_values.csv
pipeline_validation_report.txt
```

Il faut ensuite contrôler visuellement :

- le recalage CT–Tmax ;
- la normalisation ;
- l’absence d’inversion gauche–droite ;
- la couverture des deux thalamus ;
- les masques de validité.

## 7. Comparer aux résultats historiques

Avant publication du code :

- vérifiez que les 62 patients sont présents ;
- comparez les moyennes thalamiques avec le fichier ayant servi à l’analyse ;
- vérifiez l’index d’asymétrie ;
- refaites les figures ;
- documentez toute différence.

## 8. GitHub

Le fichier suivant ne doit jamais être envoyé sur GitHub :

```text
config/local_config.m
```

Le fichier `.gitignore` empêche normalement son ajout, ainsi que celui des fichiers d’imagerie et des résultats patients.

## 9. Limite importante

Les scripts ont été refactorisés à partir du code original, mais n’ont pas pu être exécutés ici avec MATLAB, SPM, CTseg et FieldTrip. Une validation locale complète est indispensable avant de créer la release GitHub `v1.0.0`.
