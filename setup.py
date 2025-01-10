from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="ivablib",
    version="0.1.0",
    author="Amelia Hercock",
    author_email="ghhercock@gmail.com",
    description="A tool for analyzing drug-induced Torsades de Pointes (TdP) risk",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Amelia3141/TdPRiskPredictor",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'ivablib': ['*.py', '*.json', '*.txt', 'data/*.csv'],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
    ],
    python_requires=">=3.9",
    install_requires=[
        "streamlit>=1.29.0",
        "plotly>=5.18.0",
        "pandas>=1.5.3",
        "biopython>=1.81",
        "requests>=2.28.2",
        "beautifulsoup4>=4.12.2",
    ],
)
