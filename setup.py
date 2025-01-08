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
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
    ],
    python_requires=">=3.10",
    install_requires=[
        "numpy>=1.24.3",
        "pandas>=2.0.3",
        "torch>=2.0.1",
        "transformers>=4.31.0",
        "sentence-transformers>=2.2.2",
        "beautifulsoup4>=4.12.2",
        "requests[socks]>=2.31.0",
        "PySocks>=1.7.1",
        "PyPDF2>=3.0.1",
        "pdfplumber>=0.10.2",
        "accelerate>=0.21.0",
        "einops>=0.6.1",
        "biopython>=1.81",
        "huggingface-hub>=0.14.1",
        "arxiv>=2.0.0",
        "plotly>=5.18.0",
        "dash>=2.14.2",
        "selenium>=4.15.2",
        "webdriver_manager>=4.0.1",
        "unpywall>=0.2.3",
        "flask>=2.3.3",
        "flask-cors>=4.0.0",
        "gunicorn>=21.2.0",
        "python-dotenv>=1.0.0"
    ],
    include_package_data=True,
    package_data={
        "ivablib": ["data/*.csv"],
    },
)
