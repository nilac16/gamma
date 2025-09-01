from setuptools.command.build_ext import build_ext
from setuptools import setup
import numpy


class BuildExt(build_ext):
    def build_extensions(self):
        if self.compiler.compiler_type == "msvc":
            for ext in self.extensions:
                ext.extra_compile_args += ["/std:c11", "/O2"]
                ext.include_dirs += [numpy.get_include()]
        build_ext.build_extensions(self)


setup(cmdclass={"build_ext": BuildExt})
