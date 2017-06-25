import h5py
import numpy as np
import time

start = time.time()
#filepath = '/phys/groups/tev/scratch1/users/jk232/'
f = h5py.File("../"+'gjj_Variables.hdf5', 'r')

'''#for testing
f = h5py.File(filepath+'mygjj_Variables.hdf5', 'r')'''

high = f['high_input'][0:10000000]
y = f['y_input'][0:10000000]

#masking nan
print 'masking nan...'

#for high
mask = ~np.isnan(high).any(axis=2)
high_input = high[mask[:,0],...]
y_input = y[mask[:,0],...]

'''for medium
maskk = ~np.isnan(high).any(axis=2)
maskk2 = ~np.any(~maskk, axis=1)
??? is this right?
'''

# Normalization
print 'Normalizing...'

for i in range(high.shape[2]):
    var = high_input[:,:,i]
    var = (var - np.mean(var))/np.std(var)
    high_input[:,:,i] = var

# Train Test split and random shuffle
print 'Train Test Split and shuffle...'

n_samples = high_input.shape[0]
tt_split = 0.8

shuff_index = np.arange(n_samples)
np.random.shuffle(shuff_index)

Train_Range = shuff_index[0:int(tt_split * n_samples)]
Test_Range = shuff_index[int(tt_split * n_samples) :]

X_train = high_input[Train_Range, :, :]
X_test = high_input[Test_Range, :, :]

y_train = y_input[Train_Range, :]
y_test = y_input[Test_Range, :]

end = time.time()
print 'prep time', (end - start)

#### Training...
print 'Building Keras models...'
from keras.layers import GRU, Highway, Dense, Dropout, MaxoutDense, Activation, Masking
from keras.models import Sequential
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras.layers.convolutional import Convolution1D, MaxPooling1D
from keras.models import Sequential
from keras.legacy.models import Graph

model = Sequential()
#model.add(Dense(25, input_dim=16))
#model.add(Masking(mask_value=-999, input_shape=(40, 2)))
model.add(GRU(25, input_shape=(1,16), dropout_W = 0.05))
# remove Maxout for tensorflow
model.add(MaxoutDense(64, 5))  #, input_shape=graph.nodes['dropout'].output_shape[1:]))
#model.add(Dense(64, activation='relu'))

model.add(Dropout(0.4))

#model.add(Highway(activation = 'relu'))

#model.add(Dropout(0.3))
model.add(Dense(3))

model.add(Activation('softmax'))

print('Compiling model...')

#adam = Adam(lr=1e-4)
model.compile(optimizer = 'Adam', loss = 'categorical_crossentropy', metrics=['accuracy'])
model.summary()

print ('Training:')
try:
    history = model.fit(X_train, y_train, batch_size=128,
        callbacks = [
            EarlyStopping(verbose=True, patience=20, monitor='val_loss'),
            ModelCheckpoint('-progress', monitor='val_loss', verbose=True, save_best_only=True)
        ],
    nb_epoch=5,
    validation_split = 0.2,
    show_accuracy=True
    )

except KeyboardInterrupt:
    print('Training ended early.')

print("history keys", history.history.keys())

y_hat = model.predict(X_test, batch_size=128)
print("y_hat.shape", y_hat.shape)
print("y_hat", y_hat[0:10])


'''ROC Curve (Recall)
True Positive Rate (tpr) = TP / P = TP / (TP + FN)
False Positive Rate (fpr) = FP / N = FP / (FP + TN) = 1 - TN / (FP + TN)

for classification:
TPR = P(test positive | is bottom jet)
FPR = P(test negative | not bottom jet)

Finding TPR:
Set probability Threshold
collect all the true signal in y
for corresponding y_hat, see how many y_hats are above Threshold'''

tpr = []
fpr = []

for i in range(100):
    th = i/float(100)
    TP = np.sum((y_hat[:,2] >= th) * y_test[:,2])
    tpr.append( TP / float(np.sum(y_test[:,2])) )

    TN = np.sum((y_hat[:,2] < th) * (1-y_test[:,2]))
    fpr.append( 1 - TN / float(np.sum(y_test[:,0]+y_test[:,1])) )

tpr = np.concatenate([[0.0], tpr])
fpr = np.concatenate([[0.0], fpr])

tprc = []
fprc = []

for i in range(100):
    th = i/float(100)
    TP = np.sum((y_hat[:,1] >= th) * y_test[:,1])
    tprc.append( TP / float(np.sum(y_test[:,1])) )

    TN = np.sum((y_hat[:,1] < th) * (1-y_test[:,1]))
    fprc.append( 1 - TN / float(np.sum(y_test[:,0]+y_test[:,2])) )

tprc = np.concatenate([[0.0], tprc])
fprc = np.concatenate([[0.0], fprc])

np.savetxt("tpr.csv", np.sort(tpr), delimiter=',')
np.savetxt("fpr.csv", np.sort(fpr), delimiter=',')
np.savetxt("tprc.csv", np.sort(tprc), delimiter=',')
np.savetxt("fprc.csv", np.sort(fprc), delimiter=',')
