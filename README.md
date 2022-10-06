# ML_proteases_cleave_Sprotein
This program creates a trained ML model for four proteases (HAT, DESC1, TMPRSS2, TMPRRSS4) known to cleave the spike protein. 

I determined parameters such as the kernel size, the number of convolutional layers, the number of neurons in the hidden layer, the batch size, and the Epochs by depth-first search in the full connection model and the CNN model. As an index for determining the parameters, I adopted the highest prediction value (accuracy) for each parameter by arranging the verification data.  

The degree of cleavage was examined by applying the ML model created with the determined parameters to the spike protein cleavage site.
