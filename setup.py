from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="tdp_risk_predictor",
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
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    include_package_data=True,
    package_data={
        "tdp_risk_predictor": ["data/*.csv"],
    },
)
