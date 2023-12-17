import keras
import PIL

# segregating the data ---------------------------------------------------------------------------------------------------------------
import os, shutil

#Path to the directory where the original dataset was uncompressed
original_dataset_dir = "C:/Users/Qba Liu/Documents/STUDIA/BIOINF_MASTER_BERLIN/SEMESTER_I/INTRODUCTION_TO_FOCUS_AREAS/DATA_SCIENCE/LAB/BreaKHis_orig"

#Directory where you’ll store your smaller dataset
# only run it once
base_dir = "C:/Users/Qba Liu/Documents/STUDIA/BIOINF_MASTER_BERLIN/SEMESTER_I/INTRODUCTION_TO_FOCUS_AREAS/DATA_SCIENCE/LAB/BreaKHis_small"
#os.mkdir(base_dir)

#Directories for the training, validation, and test splits
# only run this once
train_dir = os.path.join(base_dir, 'train') 
#os.mkdir(train_dir)
validation_dir = os.path.join(base_dir, 'validation') 
#os.mkdir(validation_dir)
test_dir = os.path.join(base_dir, 'test') 
#os.mkdir(test_dir)

#Directory with training cat pictures and dogs picture (our cats and dogs are benign and malignant)
# only run this once
train_bening_dir = os.path.join(train_dir, 'bening') 
#os.mkdir(train_bening_dir) 
train_malignant_dir = os.path.join(train_dir, 'malignant') 
#os.mkdir(train_malignant_dir)

#Directory with validation cat pictures and dog pictures (benign and malignant)
# only run it once
validation_bening_dir = os.path.join(validation_dir, 'bening') 
#os.mkdir(validation_bening_dir)
validation_malignant_dir = os.path.join(validation_dir, 'malignant') 
#os.mkdir(validation_malignant_dir)

#Directory with test cat pictures and dog pictures
# only run it once
test_benign_dir = os.path.join(test_dir, 'benign') 
#os.mkdir(test_benign_dir) 
test_malignant_dir = os.path.join(test_dir, 'malignant') 
#os.mkdir(test_malignant_dir)


#Copies the first 800 benign images to train_benign_dir
# only run this once
'''
file_names = ['benign ({}).png'.format(i) for i in range(1,800)] 
for file_name in file_names: 
 src = os.path.join(original_dataset_dir, file_name) 
 dst = os.path.join(train_bening_dir, file_name) 
 shutil.copyfile(src, dst)
'''

#Copies the next 300 benign images to validation_benign_dir
# only run once
'''
file_names = ['benign ({}).png'.format(i) for i in range(800, 1100)]
for file_name in file_names: 
 src = os.path.join(original_dataset_dir, file_name) 
 dst = os.path.join(validation_bening_dir, file_name) 
 shutil.copyfile(src, dst)
'''


#Copies the next 300 benign images to test_benign_dir
# only run this once
'''
file_names = ['benign ({}).png'.format(i) for i in range(1100, 1400)] 
for file_name in file_names: 
 src = os.path.join(original_dataset_dir, file_name) 
 dst = os.path.join(test_benign_dir, file_name) 
 shutil.copyfile(src, dst)
'''


#Copies the first 3000 malignant images to train_malignant_dir
# only run this once
"""
file_names = ['malignant ({}).png'.format(i) for i in range(1,3000)] 
for file_name in file_names: 
 src = os.path.join(original_dataset_dir, file_name) 
 dst = os.path.join(train_malignant_dir, file_name) 
 shutil.copyfile(src, dst)
"""


#Copies the next 1200 malignant images to validation_malignant_dir
# only run once
"""
file_names = ['malignant ({}).png'.format(i) for i in range(3000, 4200)] 
for file_name in file_names: 
 src = os.path.join(original_dataset_dir, file_name) 
 dst = os.path.join(validation_malignant_dir, file_name) 
 shutil.copyfile(src, dst)
"""
#Copies the next 1200 malignant images to test_malignant_dir
# only run this once
"""
file_names = ['malignant ({}).png'.format(i) for i in range(4200, 5400)] 
for file_name in file_names: 
 src = os.path.join(original_dataset_dir, file_name) 
 dst = os.path.join(test_malignant_dir, file_name) 
 shutil.copyfile(src, dst)
"""

# -------------------------------------------------------------------------------------------------------------------------------------

# building the model ---------------------------------------------------------------------------------------------------------------

from keras import layers
from keras import models
model1 = models.Sequential()
model1.add(layers.Conv2D(32, (3, 3), activation='relu',
 input_shape=(150, 150, 3)))
model1.add(layers.MaxPooling2D((2, 2)))
model1.add(layers.Conv2D(64, (3, 3), activation='relu'))
model1.add(layers.MaxPooling2D((2, 2)))
model1.add(layers.Conv2D(128, (3, 3), activation='relu'))
model1.add(layers.MaxPooling2D((2, 2)))
model1.add(layers.Conv2D(128, (3, 3), activation='relu'))
model1.add(layers.MaxPooling2D((2, 2)))
model1.add(layers.Flatten())
model1.add(layers.Dense(512, activation='relu'))
model1.add(layers.Dense(1, activation='sigmoid'))

# --------------------------------------------------------------------------------------------------------------------------------------

# optimizing the model -----------------------------------------------------------------------------------------------------------------

from keras import optimizers
model1.compile(loss='binary_crossentropy',
 optimizer=optimizers.RMSprop(learning_rate=1e-4),
 metrics=['acc'])

# --------------------------------------------------------------------------------------------------------------------------------------

# input data (image) preprocessing ----------------------------------------------------------------------------------------------------

from keras.preprocessing.image import ImageDataGenerator

#Rescales all images by 1/255
train_datagen = ImageDataGenerator(rescale=1./255)
test_datagen = ImageDataGenerator(rescale=1./255)


train_generator = train_datagen.flow_from_directory(
 train_dir,                                             #Target directory (training dataset in this case)
 target_size=(150, 150),                                #Resizes all images to 150 × 150 to make RGB to black and white
 #batch_size=20,                                        
 class_mode='binary')                                 #Because you use binary_crossentropy loss, you need binary labels.
 
validation_generator = test_datagen.flow_from_directory(
 validation_dir,
 target_size=(150, 150),
 #batch_size=20,
 class_mode='binary')

# ------------------------------------------------------------------------------------------------------------------------------------

# training and validating the model ----------------------------------------------------------------------------------------------------

history = model1.fit(
 train_generator,
 #steps_per_epoch = 5,   # this is not declared, to make keras define it as it wants to 
 epochs=5,                 
 validation_data=validation_generator,
 #validation_steps = 100,
 verbose = 2)

# ------------------------------------------------------------------------------------------------------------------------------------


# plotting the performance at each step ---------------------------------------------------------------------------------------------

from matplotlib import pyplot as plt
acc = history.history['acc']
val_acc = history.history['val_acc']
loss = history.history['loss']
val_loss = history.history['val_loss']
epochs = range(1, len(acc) + 1)

import matplotlib.pyplot as plt
acc = history.history['acc']
val_acc = history.history['val_acc']
loss = history.history['loss']
val_loss = history.history['val_loss']
epochs = range(1, len(acc) + 1)

plt.plot(epochs, acc, 'bo', label='Training acc')
plt.plot(epochs, val_acc, 'b', label='Validation acc')
plt.title('Training and validation accuracy')
plt.legend()
plt.figure()
plt.plot(epochs, loss, 'bo', label='Training loss')
plt.plot(epochs, val_loss, 'b', label='Validation loss')
plt.title('Training and validation loss')
plt.legend()
plt.show()

"""
plt.plot(epochs, acc, 'bo', label='validation acc')
plt.title('Accuracy')
plt.legend()
plt.figure()
plt.plot(epochs, loss, 'bo', label='Validation loss')
plt.title('Validation loss')
plt.legend()
plt.show()

"""
# ----------------------------------------------------------------------------------------------------------------------------------
