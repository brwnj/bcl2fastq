from setuptools import setup


setup(
    name='bcl2fastq',
    version='1.3.0',
    url='http://github.com/brwnj/bcl2fastq',
    license='MIT',
    author='Joe Brown',
    author_email='brwnjm@gmail.com',
    description='NextSeq 1k/2k specific bcl2fastq wrapper.',
    long_description=__doc__,
    py_modules=['bcl2fastq'],
    install_requires=[
        'click>=4.0',
        'pandas',
        'seaborn',
        'matplotlib',
    ],
    entry_points='''
        [console_scripts]
        bcl_to_fastq=bcl2fastq:bcl2fastq
    '''
)
