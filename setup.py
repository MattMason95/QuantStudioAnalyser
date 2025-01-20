from setuptools import setup, find_packages

setup(
    name="QuantStudioAnalyser",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        'pandas>=',
        'numpy>=',
        'matplotlib>=',
        'seaborn>=' 
    ],
    author="Matthew Mason",
    author_email="m.mason95@outlook.com",
    description="An importable library for the automated analysis of qPCR files output from QuantStudio software.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/MattMason95/",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)