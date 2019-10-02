# rtfmri-project
A Python package for: 
- taking rtfMRI data that has been preprocessed into a blood-oxygen level correlation matrix between brain regions,
- turning those correlation matrices into different types of networks, and
- measuring network properties, most of which have interpretations for the brain.

This packages focuses specifically on connectivity in the CEN (Central Executive Network) and the DMN (Default Mode Network), but it doesn't have to. The original brain scan data was taken from a study in which schizophrenic subjects (with a neurotypical control) were taught mindfulness meditation techniques, which have been proven to lower DMN connectivity and increase CEN connectivity. By using rtfMRI (real time fMRI) scans, subjects were able to play a game where they moved a dot on a screen, where the dot measured the relative connectivity between the two networks. A more thorough discussion can be found at [the study's website](https://cos.northeastern.edu/whitfield-gabrieli/projects/real-time-fmri-neurofeedback/).
