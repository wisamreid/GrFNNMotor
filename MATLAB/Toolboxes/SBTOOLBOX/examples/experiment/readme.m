% Example for the experiment description

% load the model (have a look at the model.txt file
model = SBmodel('model.txt')

% load the experiment
experiment = SBexp('experiment.exp')

% merge model and experiment
expmodel = SBmergemodexp(model,experiment)

% Have a look at the merged model
SBedit(expmodel)