# group individual replicas into sequences of 3-replicas
# Zilin Song, 23 AUG 2021
# 

import numpy, iomisc, sys
from tensorflow import keras
from tensorflow.keras import layers, optimizers, regularizers, initializers

# -------------------------
# Preprocessing.
# -------------------------
ds_idx = int(sys.argv[1])       # which dataset.
outdir = str(sys.argv[2])       # output directory

# sanity check
if ds_idx not in [0, 1, 2, 3, 4, 5]: raise ValueError('Bad ds_idx spec.\n')

sysnames  = ['toho_amp', 'toho_cex', 'toho_amp', 'toho_cex', 'both', 'both']
pathnames = ['r1ae',     'r1ae',     'r2ae',     'r2ae',     'r1ae', 'r2ae']

# Load ds
x, xpoh, y = iomisc.load_ds(sysnames[ds_idx], pathnames[ds_idx])

print('   x dim:', x.shape)
print('xpoh dim:', xpoh.shape)
print('   y dim:', y.shape)

# -------------------------
# Model definition.
# -------------------------
l_inpt_x = layers.Input(shape=(x.shape[1]))

l_dens_1 = layers.Dense(256, activation='relu')(l_inpt_x)
l_dens_1 = layers.Dropout(0.1)(l_dens_1)

l_dens_2 = layers.Dense(256, activation='relu')(l_dens_1)
l_dens_2 = layers.Dropout(0.1)(l_dens_2)

l_dens_3 = layers.Dense(256, activation='relu')(l_dens_2)
l_dens_3 = layers.Dropout(0.1)(l_dens_3)

l_inpt_xpoh = layers.Input(shape=(xpoh.shape[1]))
l_concat = layers.Concatenate()([l_dens_3, l_inpt_xpoh])

l_dens_fin = layers.Dense(256+xpoh.shape[1], activation='relu')(l_concat)
l_dens_fin = layers.Dropout(0.1)(l_dens_fin)

l_outpt    = layers.Dense(1)(l_dens_fin)

# -------------------------
# Model compile.
# -------------------------

model = keras.Model(inputs=[l_inpt_x, l_inpt_xpoh], outputs=[l_outpt])

model.compile(
    optimizer='adam', 
    loss='mse', 
    metrics=['mse', 'mae']
)

model.summary()

keras.utils.plot_model(model, to_file=f'{outdir}/model_architect.png', show_shapes=True, show_layer_names=True, dpi=300)

# -------------------------
# Model training.
# -------------------------

# Callbacks during training.
cb_checkpoint = keras.callbacks.ModelCheckpoint(filepath='./{0}/model.chkpt_tf'.format(outdir, ), )
cb_earlystopping = keras.callbacks.EarlyStopping(monitor='loss', patience=300, restore_best_weights=True)

# fitting.
history = model.fit(x=(x, xpoh), 
                    y=(y), 
                    batch_size=25,
                    epochs=300, 
                    callbacks=[cb_checkpoint, cb_earlystopping, ],
                    validation_split=0.,          
                    shuffle=True,                  # shuffle='batch' produces shitty results: don't use;
                    initial_epoch=0,               # previous epoch number.
)

# -------------------------
# Metrics
# -------------------------

y_pred = model.predict((x, xpoh))

results = model.evaluate((x, xpoh), y)

print("loss, mse, mae:", results)

model.save(f'{outdir}/model.h5')
numpy.save(f'{outdir}/y_pred.npy', y_pred)

y      = y.flatten()
y_pred = y_pred.flatten()

import math

mse = numpy.square(y - y_pred).mean()
rmse = math.sqrt(mse)
print('mse: ', mse, '\nrmse: ', rmse)

mae = numpy.abs(y-y_pred).mean()
print('mae: ', mae)
