from setuptools import setup


setup(name='TEXTURElib',
      version='0.0',
      description='Collection of texture related modules',
      author='Youngung Jeong',
      author_email='youngung.jeong@gmail.com',
      packages=['TX'],
      package_dir={'TX':'src'},
      install_requires=['fortranformat','pandas','matplotlib','numpy']
      )
      # cmdclass=cmdclass,
      #ext_modules=ext_modules)



# setup(name='TEXTURElib',
#       version='0.0',
#       description='Collection of texture related modules',
#       author='Youngung Jeong',
#       author_email='youngung.jeong@gmail.com',
#       packages=['TX'],
#       package_dir={'TX':'src'},
#       install_requires=['fortranformat','pandas','matplotlib','numpy'],
#       cmdclass=cmdclass,
#       ext_modules=ext_modules)

# ## Fortran subroutines with f2py
# ext_modules = []
# ext_modules += [
#     Extension(
#         name="pf_for_lib",
#         sources=['src/for_lib.for','src/cs.f'],
#         ## extra_compile_args=['-O3']
#         # f2py_options=['--fcompiler=intel']
#     )]
