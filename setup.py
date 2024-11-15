from setuptools import setup, find_packages

setup(
    name="methytools",
    packages=find_packages(),
    version="0.0.6",
    author="Department of research and development, Zhejiang Gaomei Genomics",
    author_email="it@gomicsgene.com",
    description="methytools is a BS-seq analysis tool suite",
    long_description="",
    long_description_content_type="text/markdown",
    license='MIT',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX :: Linux',
    ],
    package_data={
        'methytools.ref.hg38': ['*.bed.gz'],
    },
    install_requires=['pandas', 'pysam'],
    entry_points={
        'console_scripts': [
            'methytools=methytools.methytools:main',
            'mcomppost=methytools.mcomppost:main',
            'homerhelper=methytools.homerhelper:main',
        ],
    },
    py_modules=[],
    python_requires='>=3.10, <4',
)
