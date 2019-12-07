
from cffi import FFI

def compile_sniclib():
	ffibuilder = FFI()

	ffibuilder.cdef(open("snic.h").read())

	ffibuilder.set_source(
		"_snic",
		r"""
			#include "snic.h"
		""",
		sources=["snic.c"],
		library_dirs=['.'],
		extra_compile_args=['-O3', '-march=native', '-ffast-math'])

	ffibuilder.compile(verbose=True)



compile_sniclib()
