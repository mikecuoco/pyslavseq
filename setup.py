from setuptools import find_packages
from setuptools import setup

setup(name='pyslavseq',
      version='0.1',
      description='SLAV-seq analysis package',
      url='http://github.com/apuapaquola/pyslavseq',
      author='Apuã Paquola',
      author_email='apuapaquola@gmail.com',
      license='GPLv3',
      packages=find_packages('src'),
      package_dir={'': 'src'},
      include_package_data=True,
      entry_points={
          'console_scripts': [
              'pyslavseq_extract_features = pyslavseq.features.get_window_features_occupied:main'
          ]
      },
      zip_safe=False)
