from setuptools import setup, find_packages

setup(
    name="align",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "typing"
    ],
    python_requires=">=3.7",
)