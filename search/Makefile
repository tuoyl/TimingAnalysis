target:
	g++ FreqSearch.cxx readparfile.cxx /Users/tuoyouli/cpplibrary/MyFun.cxx -I${HEADAS}/include -L${HEADAS}/lib -std=c++11 -I./include -I/Library/Python/2.7/site-packages/numpy/core/include/ -I/usr/local/gcc_system_root/usr/include/python2.7/ -lpython  -lcfitsio  -o FreqSearch
	./FreqSearch -f ../testdata/P021405900101_HE_screen.fits -p ../testdata/RXJ0209.6.par
