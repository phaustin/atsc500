from setuptools import setup, find_packages

setup(
    name = "a500",
    packages=find_packages(),
    entry_points={
          'console_scripts': [
              'killjobs = a500.utils.killjobs:killjobs',
              'pyncdump = a500.utils.ncdump:main'
          ]
    },

)
