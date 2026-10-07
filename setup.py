from setuptools import setup, find_packages

setup(
    name='three-body-problem',
    version='0.1.0',
    description='A research-grade Python package for simulating the gravitational three-body problem',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Shreeharsh Shinde',
    author_email='shreeharsh157@gmail.com',
    license='MIT',
    url='https://github.com/shreeharshshinde/three-body-problem',
    packages=find_packages(),
    install_requires=[
        'numpy>=1.24',
        'scipy>=1.10',
        'matplotlib>=3.7',
    ],
    extras_require={
        'dev': [
            'pytest>=7.0.0',
            'pytest-cov>=4.0',
        ],
        'notebooks': [
            'jupyterlab>=4.0',
        ],
    },
    python_requires='>=3.10',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
)
