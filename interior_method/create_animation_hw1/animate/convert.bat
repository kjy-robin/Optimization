echo on

cd png

FOR %%X in ("*.png") DO ..\convert\convert.exe -resize 500x500 %%X ..\gif\%%X.gif

cd ..
