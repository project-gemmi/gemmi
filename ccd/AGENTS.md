This is gemmi - a project in C++11
I'm workin on "gemmi drg" ( ./prog/dlg.cpp )  command which should reproduce part of the CCP4 Acedrg
For reference, Acedrg code is in ./acedrg/. What we implement, is reading Acedrg tables (./acedrg/tables/) and calculating bond and angle restraints for a given CCD monomer.
We try to exactly reproduce Acedrg algos: protonation /deprotonationlogic, Acedrg atom types, etc.
Protonation, atom types and bond restraints work fine. Nos we need to sort out angles.

The actual code is in `./include/gemmi/ace*.hpp` and `./src/ace*.cpp`.
I have a test rig  for comparing out results against Acedrg, I'll be copy pasting differences for you.
Your goal is to eliminate these differences.
Data (monomers for testing) is in (using ABC as an example):
input from CCD: ccd/orig/a/ABC.cif
output from Acedrg: ./ccd/acedrg/a/ABC.cif
logs from Acedrg: ./ccd/acedrg/a/ABC_TMP/* .log
output from gemmi: ./ccd/gemmi/a/ABC.cif

To rebuild gemmi, run ./run-tests.sh
or directly:
cmake --build /home/wojdyr/gemmi/gemmi/build -j
to run gemmi drg you can use, e.g.:
 ACEDRG_TABLES=./acedrg/tables ./build/gemmi drg ccd/orig/ALA.cif /tmp/ala_gemmi.cif
then, to compare it against acedrg:
 ./build/gemmi mondiff ccd/acedrg/ALA.cif /tmp/ala_gemmi.cif

If you think that AceDRG is wrong, write a short description that can be sent to Acedrg developers.

If not sure what exact rules Acedrg is using, don't guess!!!, write a self-contained question for
Acedrg experts. when asking question, don't list answer and expect free-text answer.

The latest acedrg is installed in pythonic virtual env:
/home/wojdyr/app/acedrg/venv
but it's not reliable, so the one from CCP4 is better for our needs.
