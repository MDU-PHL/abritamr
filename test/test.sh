#!/usr/bin/env bats

setup () {
  name="abritamr"
  bats_require_minimum_version 1.5.0
  dir=$(dirname "$BATS_TEST_FILENAME")
  cd "$dir"
  exe="$name"
}

@test "Version" {
  run -0 $exe --version
  [[ "$output" =~ "$name " ]]
}
@test "Help" {
  run -0 $exe --help
  [[ "$output" =~ "detection" ]]
}
@test "No parameters" {
  run -0 $exe
  [[ "$output" =~ "usage:" ]]
}
@test "Bad option" {
  run ! $exe --doesnotexist
}
@test "Non-existent input file" {
  run ! $exe run -c doesnotexist.fasta
}
@test "Empry input file" {
  run ! $exe run -c /dev/null
}
@test "Small input file" {
  run -0 $exe run -c small.fna
  [[ -r "abritamr/amrfinder.out" ]]
}
