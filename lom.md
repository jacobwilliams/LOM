project: lom
project_dir: ./
src_dir: ./app
output_dir: ./doc
media_dir: ./media
project_github: https://github.com/jacobwilliams/LOM
summary: Low Lunar Orbit Maintenance
author: Jacob Williams
github: https://github.com/jacobwilliams
predocmark_alt: >
predocmark: <
docmark_alt:
docmark: !
display: public
         private
         protected
source: true
graph: true
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html
            pyplot_module:https://jacobwilliams.github.io/pyplot-fortran
            ddeabm_module:https://jacobwilliams.github.io/ddeabm/
            json_module:https://jacobwilliams.github.io/json-fortran/
            argv_module:https://jacobwilliams.github.io/argv-fortran/
            fortran_astrodynamics_toolkit:https://jacobwilliams.github.io/Fortran-Astrodynamics-Toolkit/

{!README.md!}