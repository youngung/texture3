from setuptools import setup
setup(name='TEXTURElib',
      version='1.0',
      description='Collection of texture related modules',
      author='Youngung Jeong',
      author_email='yjeong@changwon.ac.kr',
      packages=['TX'],
      package_dir={'TX':'src'},
      install_requires=['fortranformat','pandas','matplotlib','numpy','scipy','shapely']
      )
