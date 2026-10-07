# Li Sigma5 (310)[001] 2x2 -- Windows local LAMMPS runner
# Requires lmp.exe with OpenKIM support in PATH.
$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $PSScriptRoot
$Kim = "SNAP_ZuoChenLi_2019_Li__MO_732106099012_000"
$Seed = "12345"

Push-Location $Root
try {
  & lmp -in lammps/in.0_check -var outdir stage0_build_validate -var kim_model $Kim
  & lmp -in lammps/in.1_minimize -var outdir stage1_relax -var input_data stage0_build_validate/Sigma5_310_K6_J6_L20.data -var kim_model $Kim
  & lmp -in lammps/in.2_highT_npt -var outdir stage2_highT_NPT -var seed $Seed -var kim_model $Kim
  & lmp -in lammps/in.3_quench -var outdir stage3_quench -var kim_model $Kim
  & lmp -in lammps/in.4_10K_hold -var outdir stage4_10K_hold -var kim_model $Kim
}
finally {
  Pop-Location
}
