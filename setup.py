from setuptools import setup, find_packages

setup(name='Superharm',
      version='1.0',
      packages = find_packages(),
      description='Dynamic System Simulator',
      license='GPL2',
      author='Horacio Vargas Guzman',
      author_email='horacio.v.g@gmail.com',
      url='https://github.com/horatz/superharm',
      scripts = ['Superharm'],
      include_package_data = True,
      package_data = {
        'superharm': ['*',],
      }
)
