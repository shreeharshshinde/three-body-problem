from setuptools import setup, find_packages

setup(
    name='three-body-problem',
    version='0.1.0',
    description='A research-grade Python package for simulating the gravitational three-body problem',
    author='Your Name',
    author_email='your.email@example.com',
    license='MIT',
    packages=find_packages(),
    install_requires=[
        'numpy>=1.22.0',
        'scipy>=1.8.0',
        'matplotlib>=3.5.0',
        'jupyter>=1.0.0',
    ],
    extras_require={
        'dev': ['pytest>=7.0.0'],
    },
    python_requires='>=3.10',
)
