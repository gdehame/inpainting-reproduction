# Image inpainting

## Structure
This project folder is organized as the CGDI practicals.

The `CGDI_Proj.pdf` file contains the report.

The `/data` folder contains the test files, 
`gen-tests.sh` is a script used to run all tests.

The results of my experiments are in `/build/Proj`, 
these are named `testX.jpg` with `X` the test number.

## Build
To build the project you can open `/build` in a terminal and execute `make`.

This will generate in `/build/Proj` the executable `inpainting`.

## Run
Then you can execute the inpainter with `inpainting -s source_image -m mask_image -o output_image`.

You can also add the option `--sbs` to ouput the image at each stage (was used for debugging).
