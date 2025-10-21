from setuptools import setup, find_packages
from pathlib import Path

root = Path(__file__).parent
top_level_shims = [
    p.stem for p in root.glob('*.py')
    if p.name not in {'setup.py', '__init__.py'} and not p.name.startswith('.')
]

print('')
print('top_level_shims:')
print(top_level_shims)
print('')

setup(
    name='pythonmodules',
    version='0.1.0',
    description='Utilities',
    python_requires='>=3.8',
    packages=find_packages(include=['pythonmodules', 'pythonmodules.*']),
    py_modules=top_level_shims,  # installs the shim modules too
    include_package_data=True,
    install_requires=[
        # 'pandas>=2.0',
    ],
)

