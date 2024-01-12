#!/bin/bash

cd build
make
cd ..

./build/Proj/inpainting -s ./data/wall-robot-400.png -m ./data/wall-robot-400-mask.png -o ./build/Proj/test0.jpg
./build/Proj/inpainting -s ./data/image1.jpg -m ./data/mask1.jpg -o ./build/Proj/test1.jpg
./build/Proj/inpainting -s ./data/image2.jpg -m ./data/mask2.jpg -o ./build/Proj/test2.jpg
./build/Proj/inpainting -s ./data/image3.jpg -m ./data/mask3.jpg -o ./build/Proj/test3.jpg
./build/Proj/inpainting -s ./data/image4.jpg -m ./data/mask4.jpg -o ./build/Proj/test4.jpg
./build/Proj/inpainting -s ./data/image5.jpg -m ./data/mask5.jpg -o ./build/Proj/test5.jpg
./build/Proj/inpainting -s ./data/image6.jpg -m ./data/mask6.jpg -o ./build/Proj/test6.jpg
./build/Proj/inpainting -s ./data/image7.jpg -m ./data/mask7.jpg -o ./build/Proj/test7.jpg
