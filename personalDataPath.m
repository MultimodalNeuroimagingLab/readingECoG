function localDataPath = personalDataPath()

% function that contains local data path

localDataPath.input  = '';
localDataPath.output = '';

% path to freesurfer derivatives
localDataPath.fs     = '';

% other git repositories
addpath(genpath('/Users/xxx/Documents/git_repo/knkutils'))
addpath(genpath('/Users/xxx/Documents/git_repo/mnl_ieegBasics'))
