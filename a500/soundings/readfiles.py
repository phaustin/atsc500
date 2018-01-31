import pathlib
the_dir=pathlib.Path('soundingdir')
print(the_dir.exists())
print(list(the_dir.glob('**/*')))
