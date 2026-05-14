<#
.SYNOPSIS
  Assemble a folder tree + README for Zenodo dataset upload (large inputs not on GitHub).

.DESCRIPTION
  Copies known paths from <DataParent>\data\ into a staged directory with numbered
  subfolders and copies zenodo_dataset_bundle\README.md to the staging root.
  Optionally zips the result. See docs/ZENODO.md and zenodo_dataset_bundle\README.md.

.PARAMETER DataParent
  Directory that contains a child folder named "data" (e.g. UNDEFINED or paper-github).

.PARAMETER OutDir
  Output staging directory (will be created). Default: paper-github\zenodo_dataset_staging\

.PARAMETER Zip
  If set, creates OutDir.zip next to OutDir (removes existing zip first).

.EXAMPLE
  .\scripts\build_zenodo_dataset_bundle.ps1 -DataParent "C:\Users\you\Desktop\masters\UNDEFINED" -Zip
#>
param(
    [Parameter(Mandatory = $true)]
    [string] $DataParent,

    [Parameter(Mandatory = $false)]
    [string] $OutDir = "",

    [switch] $Zip
)

$ErrorActionPreference = "Stop"
$RepoRoot = Resolve-Path (Join-Path $PSScriptRoot "..")
if (-not $OutDir) {
    $OutDir = Join-Path $RepoRoot "zenodo_dataset_staging"
}

$DataRoot = Join-Path $DataParent "data"
if (-not (Test-Path $DataRoot)) {
    throw "Data directory not found: $DataRoot (check -DataParent points at folder containing 'data')"
}

$Netmhc = Join-Path $DataRoot "netmhc"
$ReadmeSrc = Join-Path $RepoRoot "zenodo_dataset_bundle\README.md"
if (-not (Test-Path $ReadmeSrc)) {
    throw "Missing template README: $ReadmeSrc"
}

Remove-Item -Path $OutDir -Recurse -Force -ErrorAction SilentlyContinue
New-Item -ItemType Directory -Path $OutDir -Force | Out-Null

$dirs = @(
    "01_tcga_lncrna_expression",
    "02_smprot_raw",
    "03_netmhc_merged",
    "04_netmhc_wide_xls",
    "05_iedb_inputs",
    "06_predictions_optional",
    "07_provenance",
    "08_curated_lists_optional"
)
foreach ($d in $dirs) {
    New-Item -ItemType Directory -Path (Join-Path $OutDir $d) -Force | Out-Null
}

function Copy-IfExists {
    param([string] $SourcePath, [string] $DestDir)
    if (Test-Path $SourcePath) {
        Copy-Item -LiteralPath $SourcePath -Destination $DestDir -Force
        return $true
    }
    return $false
}

$copied = [System.Collections.Generic.List[string]]::new()
$missing = [System.Collections.Generic.List[string]]::new()

# --- 01 expression ---
foreach ($name in @(
        "primary_exp_stage_lncRNA.csv",
        "primary_exp_metastasis_lncRNA.csv",
        "primary_exp_stage_final.csv",
        "primary_exp_metastasis_final.csv"
    )) {
    $p = Join-Path $DataRoot $name
    if (Copy-IfExists $p (Join-Path $OutDir "01_tcga_lncrna_expression")) {
        $copied.Add("01_tcga_lncrna_expression\$name")
    }
    else { $missing.Add($p) }
}
Get-ChildItem -Path $DataRoot -Filter "primary_exp_*_final.csv" -File -ErrorAction SilentlyContinue | ForEach-Object {
    $dest = Join-Path $OutDir "01_tcga_lncrna_expression"
    if (-not (Test-Path (Join-Path $dest $_.Name))) {
        Copy-Item $_.FullName -Destination $dest -Force
        $copied.Add("01_tcga_lncrna_expression\$($_.Name)")
    }
}

# --- 02 SmProt ---
$sm = Join-Path $DataRoot "SmProt2.txt"
if (Copy-IfExists $sm (Join-Path $OutDir "02_smprot_raw")) { $copied.Add("02_smprot_raw\SmProt2.txt") } else { $missing.Add($sm) }

# --- 03 merged TSV ---
if (Test-Path $Netmhc) {
    Get-ChildItem -Path $Netmhc -Filter "*_with_iedb.tsv" -File -ErrorAction SilentlyContinue | ForEach-Object {
        Copy-Item $_.FullName -Destination (Join-Path $OutDir "03_netmhc_merged") -Force
        $copied.Add("03_netmhc_merged\$($_.Name)")
    }
}

# --- 04 wide xls (cohort) ---
if (Test-Path $Netmhc) {
    Get-ChildItem -Path $Netmhc -Filter "netmhcpan*.xls" -File -ErrorAction SilentlyContinue | ForEach-Object {
        Copy-Item $_.FullName -Destination (Join-Path $OutDir "04_netmhc_wide_xls") -Force
        $copied.Add("04_netmhc_wide_xls\$($_.Name)")
    }
}

# --- 05 IEDB ---
if (Test-Path $Netmhc) {
    Get-ChildItem -Path $Netmhc -Filter "iedb_*.csv" -File -ErrorAction SilentlyContinue | ForEach-Object {
        Copy-Item $_.FullName -Destination (Join-Path $OutDir "05_iedb_inputs") -Force
        $copied.Add("05_iedb_inputs\$($_.Name)")
    }
    Get-ChildItem -Path $Netmhc -Filter "iedb_*.json" -File -ErrorAction SilentlyContinue | ForEach-Object {
        Copy-Item $_.FullName -Destination (Join-Path $OutDir "05_iedb_inputs") -Force
        $copied.Add("05_iedb_inputs\$($_.Name)")
    }
}

# --- 06 predictions subtree ---
$predSrc = Join-Path $Netmhc "predictions"
$predDst = Join-Path $OutDir "06_predictions_optional"
if (Test-Path $predSrc) {
    Copy-Item -Path $predSrc -Destination $predDst -Recurse -Force
    $copied.Add("06_predictions_optional\ (tree from data\netmhc\predictions)")
}

# --- 07 provenance ---
$man = Join-Path $DataRoot "MANIFEST_COPIED.txt"
if (Copy-IfExists $man (Join-Path $OutDir "07_provenance")) { $copied.Add("07_provenance\MANIFEST_COPIED.txt") } else { $missing.Add($man) }

# --- 08 optional curated ---
$geneSmall = Join-Path $DataRoot "lncrna_genes_small.csv"
if (Copy-IfExists $geneSmall (Join-Path $OutDir "08_curated_lists_optional")) {
    $copied.Add("08_curated_lists_optional\lncrna_genes_small.csv")
}

Copy-Item -LiteralPath $ReadmeSrc -Destination (Join-Path $OutDir "README.md") -Force

# Manifest
$manifestPath = Join-Path $OutDir "FILE_MANIFEST.txt"
$lines = @(
    "Generated: $(Get-Date -Format o)",
    "Source data root: $DataRoot",
    "",
    "=== Copied into bundle (relative paths) ==="
)
foreach ($c in $copied) { $lines += $c }
$lines += ""
$lines += "=== Missing at source (not an error if you omit optional files) ==="
foreach ($m in $missing) { $lines += $m }
$lines | Set-Content -Path $manifestPath -Encoding utf8

Write-Host "Staging complete: $OutDir"
Write-Host "Files/folders recorded: $($copied.Count); missing paths logged: $($missing.Count)"
if ($missing.Count -gt 0) {
    Write-Host "(Review FILE_MANIFEST.txt for missing list.)" -ForegroundColor DarkYellow
}

if ($Zip) {
    $zipPath = "$OutDir.zip"
    if (Test-Path $zipPath) { Remove-Item $zipPath -Force }
    # Compress folder contents so unzip creates numbered folders + README at top level
    Compress-Archive -Path (Join-Path $OutDir "*") -DestinationPath $zipPath -Force
    Write-Host "Zip: $zipPath"
}
