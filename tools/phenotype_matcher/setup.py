"""
Setup script for the Phenotype Matcher tool.
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="phenotype-matcher",
    version="1.0.0",
    author="PAVS Team",
    description="A standalone tool for matching clinical phenotype descriptions to ontology identifiers using Graph RAG",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/your-org/pavs",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Healthcare Industry",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "pyhpo>=3.2.0",
        "pronto>=2.5.0",
        "sentence-transformers>=2.2.0",
        "scikit-learn>=1.0.0",
        "requests>=2.28.0",
        "torch>=1.13.0",
    ],
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=4.0.0",
            "black>=23.0.0",
            "flake8>=6.0.0",
            "mypy>=1.0.0",
        ],
        "progress": [
            "tqdm>=4.65.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "phenotype-matcher=phenotype_matcher.cli:main",
            "phenotype-index=phenotype_matcher.index_manager:main",
        ],
    },
)
