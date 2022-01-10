# group individual replicas into sequences of 3-replicas
# Zilin Song
# 

import numpy, iomisc, sys
from tensorflow import keras
from tensorflow.keras import layers

# -------------------------
# Preprocessing.
# -------------------------
outdir   =   str(sys.argv[1])       # output directory
ds_idx   =   int(sys.argv[2])       # which dataset.

def test(ds_idx, do_rate, n_units):

    # sanity check
    if ds_idx not in [0, 1, 2, 3, 4, 5]: raise ValueError('Bad ds_idx spec.\n')

    sysnames  = ['toho_amp', 'toho_cex', 'toho_amp', 'toho_cex', 'both', 'both']
    pathnames = ['r1ae',     'r1ae',     'r2ae',     'r2ae',     'r1ae', 'r2ae']

    # Load ds
    xtrain, xpohtrain, ytrain, xtest, xpohtest, ytest= iomisc.load_ds(sysnames[ds_idx], pathnames[ds_idx])

    print('   xtrain dim:', xtrain.shape)
    print('xpohtrain dim:', xpohtrain.shape)
    print('   ytrain dim:', ytrain.shape)

    # -------------------------
    # Model definition.
    # -------------------------
    l_inpt_x = layers.Input(shape=(xtrain.shape[1]))

    l_dens_1 = layers.Dense(n_units, activation='relu')(l_inpt_x)
    l_dens_1 = layers.Dropout(0.1)(l_dens_1)

    l_dens_2 = layers.Dense(n_units, activation='relu')(l_dens_1)
    l_dens_2 = layers.Dropout(0.1)(l_dens_2)

    l_dens_3 = layers.Dense(n_units, activation='relu')(l_dens_2)
    l_dens_3 = layers.Dropout(0.1)(l_dens_3)

    l_inpt_xpoh = layers.Input(shape=(xpohtrain.shape[1]))
    l_concat = layers.Concatenate()([l_dens_3, l_inpt_xpoh])

    l_dens_fin = layers.Dense(n_units+xpohtrain.shape[1], activation='relu')(l_concat)
    l_dens_fin = layers.Dropout(0.1)(l_dens_fin)

    l_outpt    = layers.Dense(1)(l_dens_fin)

    # -------------------------
    # Model compile.
    # -------------------------

    # model = keras.Model(inputs=[l_inpt_rcrr, l_inpt_y2], outputs=[ly_outp])
    model = keras.Model(inputs=[l_inpt_x, l_inpt_xpoh], outputs=[l_outpt])

    model.compile(
        optimizer='adam', 
        loss='mse', 
        metrics=['mse', 'mae']
    )

    model.summary()

    # -------------------------
    # Model training.
    # -------------------------

    # Callbacks during training.
    cb_checkpoint = keras.callbacks.ModelCheckpoint(filepath='./{0}/model.chkpt_tf'.format(outdir, ), )
    cb_earlystopping = keras.callbacks.EarlyStopping(monitor='loss', patience=300, restore_best_weights=True, )

    # fitting.
    history = model.fit(x=(xtrain, xpohtrain), 
                        y=(ytrain), 
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

    ypred = model.predict((xtest, xpohtest))

    results = model.evaluate((xtest, xpohtest), ytest)

    print("loss, mse, mae:", results)

    ytest = ytest.flatten()
    ypred = ypred.flatten()

    import math

    test_mse  = numpy.square(ytest - ypred).mean()
    test_rmse = math.sqrt(test_mse)
    test_mae  = numpy.abs(ytest - ypred).mean()
    
    ypred = model.predict((xtrain, xpohtrain))
    results = model.evaluate((xtrain, xpohtrain), ytrain) 
    
    ytrain = ytrain.flatten()
    ypred = ypred.flatten()
    
    train_mse  = numpy.square(ytrain - ypred).mean()
    train_rmse = math.sqrt(train_mse)
    train_mae  = numpy.abs(ytrain - ypred).mean()

    with open(f'{outdir}/performance.log', 'a') as fo:
        fo.write(f'nlay = {n_units}\tdo_rate={do_rate}\ttest_mse: {test_mse:10.4f}\ttest_rmse: {test_rmse:10.4f}\ttest_mae: {test_mae:10.4f}\t')
        fo.write(f'train_mse: {train_mse:10.4f}\ttrain_rmse: {train_rmse:10.4f}\ttrain_mae: {train_mae:10.4f}\n')
    
for do_rate in [0., 0.1, 0.2, 0.3, 0.4]:
    for n_units in [64, 128, 256, 512]:
        test(ds_idx, do_rate, n_units)
