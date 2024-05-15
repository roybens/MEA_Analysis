from setuptools import setup, find_packages

# with open("README.md", "r") as fh:
#     long_description = fh.read()

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='mea_processing_library',
    version='0.0',
    author='Adam Weiner',
    author_email='amwe@ucdavis.edu',
    description='Under Construction',
    # long_description=long_description,
    # long_description_content_type='text/markdown',
    # url='https://github.com/yourusername/mea_processing_library',
    packages=find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        # 'License :: OSI Approved :: MIT License',
        # 'Operating System :: OS Independent',
    ],
    install_requires=required,
)