import numpy as np # Library imports
import matplotlib.pyplot as plt
import pandas as pd
import keras
import os
from sklearn.preprocessing import MinMaxScaler
import tensorflow as tf
from scipy.interpolate import PchipInterpolator

# Matplotlib style
plt.style.use('dark_background')

# GPU parallelisation (requires python 3.10 and CUDA 8.0 as of right now, will run as normal on CPU if not possible)
gpus = tf.config.experimental.list_physical_devices('GPU')
if gpus:
    for gpu in gpus:
        tf.config.experimental.set_memory_growth(gpu, True)

# ---------------------------------------------------------------------------
# Data Paths - Some data files contain anomalies, these have been removed from these lists
# ---------------------------------------------------------------------------

# Cell 12 - 88.55% SoH (test data is witheld from dataset to stop data leakage)
directory_cell12 = r'CELL 12 DIRECTORY' # Change directory !!

data_train_list_cell12 = ['00.txt','05.txt','10.txt','15.txt',
    '20.txt','50.txt','60.txt','65.txt',
    '70.txt','75.txt','80.txt','85.txt']

data_test_list_cell12 = ['40.txt','55.txt']

drt_train_list_cell12 = ['00DRT.csv','05DRT.csv','10DRT.csv','15DRT.csv',
    '20DRT.csv','50DRT.csv','60DRT.csv','65DRT.csv',
    '70DRT.csv','75DRT.csv','80DRT.csv','85DRT.csv']

drt_test_list_cell12 = ['40DRT.csv','55DRT.csv']

soc_train_cell12 = [0.00,5.65,11.29,16.94,22.59,
    56.47,67.76,
    73.41,79.06,
    84.71,90.35,96.00]

soc_test_cell12 = [45.17,62.12]

soh_train_cell12 = 88.55

soh_test_cell12 = 88.55

# Cell 1 - 85.94% SoH
directory_cell1 = r'CELL 1 DIRECTORY' # Change directory !!

data_train_list_cell1 = ['00.txt','05.txt','15.txt',
    '20.txt','25.txt','40.txt',
    '60.txt','70.txt','75.txt','80.txt','85.txt']

data_test_list_cell1 = ['35.txt','55.txt']

drt_train_list_cell1 = ['00DRT.csv','05DRT.csv',
    '15DRT.csv','20DRT.csv',
    '25DRT.csv','35DRT.csv','55DRT.csv',
    '60DRT.csv',
    '70DRT.csv','75DRT.csv','80DRT.csv',
    '85DRT.csv']

drt_test_list_cell1 = ['30DRT.csv','40DRT.csv','50DRT.csv']

soc_train_cell1 = [0.00,5.82,17.45,23.27,
    29.09,40.72,64.00,69.81,
    81.45,87.27,93.09,98.91]

soc_test_cell1 = [34.91,46.54,58.18]

soh_train_cell1 = 85.94

soh_test_cell1 = 85.94

# Old cell 1 - 100% SoH
directory_cell01 = r'HEALTHY CELL DIRECTORY' # Change directory !!

data_train_list_cell01 = ['00.txt','20.txt','30.txt','40.txt','50.txt',
    '60.txt','80.txt','90.txt','100.txt']

data_test_list_cell01 = ['70.txt']

drt_train_list_cell01 = ['00DRT.csv','20DRT.csv',
    '30DRT.csv','40DRT.csv',
    '50DRT.csv','60DRT.csv',
    '80DRT.csv','90DRT.csv',
    '100DRT.csv']

drt_test_list_cell01 = ['70DRT.csv']

soc_train_cell01 = [0.00,25.4,34.8,44.1,53.4,
    62.7,81.4,90.7,100.0]

soc_test_cell01 = [72.0]

soh_train_cell01 = 100

soh_test_cell01 = 100

# Final Test Dataset - Cell 10 83.30% SoH
directory_cell10 = r'CELL 10 DIRECTORY' # Change directory !!

data_cell10 = ['00.txt','10.txt',
    '20.txt','30.txt','35.txt','40.txt',
    '45.txt','55.txt','60.txt','65.txt',
    '70.txt','85.txt']

drt_cell10 = ['00DRT.csv',
    '10DRT.csv','20DRT.csv',
    '30DRT.csv',
    '40DRT.csv','45DRT.csv',
    '55DRT.csv','60DRT.csv',
    '65DRT.csv',
    '70DRT.csv']

soc_cell10 = [0.00,12.00,24.00,
    36.01,42.01,
    48.01,54.01,66.01,72.01,
    78.01,84.01]

soh_cell10 = 83.30

# ---------------------------------------------------------------------------
# Function definitions and helpers
# ---------------------------------------------------------------------------

def construct_feature(directory, file_name, drt_file_name, soc, visuals=False):
    images = []

    # EIS data
    path = os.path.join(directory, file_name)
    data = pd.read_csv(path, sep="\t")
    freq = data.iloc[:,0].to_numpy()
    ZI = data.iloc[:, 1].to_numpy()
    ZII = data.iloc[:, 2].to_numpy()

    # Remove invalid impedance region with mask
    mask = -ZII > 0
    ZI = ZI[mask]
    ZII = ZII[mask]
    freq = freq[mask]
    
    # DRT data
    path_drt = os.path.join(directory, drt_file_name)
    drt_data = pd.read_csv(path_drt)
    freq_drt = drt_data.iloc[:, 0].to_numpy()
    gamma = drt_data.iloc[:, 1].to_numpy()
    
    tau = 1 / freq
    tau_drt = 1 / freq_drt
    tau = np.log(tau)
    tau_drt = np.log(tau_drt)
    # freq =  np.log(freq)
    freq_drt = np.log(freq_drt)

    # PCHIP interpolation - Shape preserving interpolator (piecewise cubic Hermite interpolating polynomial), 
    # extends EIS data domain to fit number of DRT points
    ZI_spline = PchipInterpolator(tau, ZI) 
    ZII_spline = PchipInterpolator(tau, ZII)
    gamma_spline = PchipInterpolator(tau_drt,gamma)
    tau_min = max(tau.min(),tau_drt.min())
    tau_max = min(tau.max(),tau_drt.max())
    tau_new = np.linspace(tau_min,tau_max,510) 
    ZI_new = ZI_spline(tau_new) 
    ZII_new = ZII_spline(tau_new)
    gamma_new = gamma_spline(tau_new)

    # Append default image (EIS + DRT)
    image = np.array([ZI_new,ZII_new,gamma_new,tau_new])
    # image = np.array([ZI_new,ZII_new]) # EIS version only
    images.append(image)

    # Data augmentation - copy original set or it will be adapted in place and create additions
    shift_ZI = np.random.uniform(-0.001,0.001,1)
    shift_ZII = np.random.uniform(-0.001,0.000,1)
    ZI_original = ZI_new.copy()
    ZII_original = ZII_new.copy()

    # Plot visualisation (setting visuals==True will show every image!)
    # if visuals:
    #     plt.figure(figsize=(10, 7))
    #     plt.plot(ZI,-ZII,label='Original')
    #     plt.plot(ZI_new,-ZII_new,label='PCHIP Interpolation')
    #     plt.xlabel(r'$Z_I$')
    #     plt.ylabel(r'$-Z_{II}$')
    #     plt.title(f"Nyquist - SOC: {soc}")
    #     plt.legend()
    #     # plt.grid(alpha=0.3)
    #     plt.show()

        # plt.figure(figsize=(10, 7))
        # plt.semilogx(10**tau_drt,gamma,label='Original')
        # plt.semilogx(10**tau_new,gamma_new,label='Fit')
        # plt.xlabel(r'$\tau$')
        # plt.ylabel(r'$\gamma(\tau)$')
        # plt.title(f"DRT - SOC: {soc}")
        # plt.legend()
        # plt.grid(alpha=0.3)
        # plt.show()

        # plt.figure(figsize=(10, 7))
        # plt.title(f"SOC = {soc}")
        # plt.plot(ZI,-ZII,'--',label='Original')
        # plt.plot(ZI_original,-ZII_original,label='PCHIP Interpolation')
        # plt.xlabel(r'$Z_I$')
        # plt.ylabel(r'$-Z_{II}$')
        # plt.legend()
        # plt.grid(alpha=0.3)

    # Data augmentation
    for a, b in zip(shift_ZI, shift_ZII):
        ZI_shifted = ZI_original + a
        ZII_shifted = ZII_original + b

        image = np.array([ZI_shifted,ZII_shifted,gamma_new,tau_new])
        # image = np.array([ZI_shifted,ZII_shifted]) # EIS version only
        images.append(image)

        # if visuals: plt.plot(ZI_shifted,-ZII_shifted)

    # if visuals: plt.show()

    # output (y)
    output = np.repeat(soc,repeats=len(images))

    return images, output

# ---------------------------------------------------------------------------
# Dataset Construction
# ---------------------------------------------------------------------------

# Create Dataset
x_train = []
x_test = []
y_train = []
y_test = []

# Apply data augmentation to each dataset
for file,drt,s in zip(data_train_list_cell12,drt_train_list_cell12,soc_train_cell12):
    images,output = construct_feature(directory_cell12,file,drt,s,soh_train_cell12)
    x_train.extend(images)
    y_train.extend(output)

for file,drt,s in zip(data_test_list_cell12,drt_test_list_cell12,soc_test_cell12):
    images,output = construct_feature(directory_cell12,file,drt,s,soh_test_cell12)
    x_test.extend(images)
    y_test.extend(output)

for file,drt,s in zip(data_train_list_cell1,drt_train_list_cell1,soc_train_cell1):
    images,output = construct_feature(directory_cell1,file,drt,s,soh_train_cell1)
    x_train.extend(images)
    y_train.extend(output)

for file,drt,s in zip(data_test_list_cell1,drt_test_list_cell1,soc_test_cell1):
    images,output = construct_feature(directory_cell1,file,drt,s,soh_test_cell1)
    x_test.extend(images)
    y_test.extend(output)

for file,drt,s in zip(data_train_list_cell01,drt_train_list_cell01,soc_train_cell01):
    images,output = construct_feature(directory_cell01,file,drt,s,soh_train_cell01)
    x_train.extend(images)
    y_train.extend(output)

for file,drt,s in zip(data_test_list_cell01,drt_test_list_cell01,soc_test_cell01):
    images,output = construct_feature(directory_cell01,file,drt,s,soh_test_cell01)
    x_test.extend(images)
    y_test.extend(output)

# Stack and transpose datasets (create meaningful images, can be adapted to produce different images)
x_train = np.stack(x_train)
x_test = np.stack(x_test)
y_train = np.array(y_train)
y_test = np.array(y_test)

x_train = np.transpose(x_train, (0, 2, 1))
x_test = np.transpose(x_test, (0, 2, 1))
X_train = x_train[..., np.newaxis]
X_test = x_test[..., np.newaxis]

# Input scaling
scalers = []

for channel in range(4): # Change to 2 if only EIS version of image
    scaler = MinMaxScaler(feature_range=(0,1))

    # Training channel (fit on training data)
    X_channel = X_train[:, :, channel, 0]
    scaler.fit(X_channel)

    # Transform
    X_train[:, :, channel, 0] = scaler.transform(X_channel)
    X_test[:, :, channel, 0] = scaler.transform(X_test[:, :, channel, 0])
    scalers.append(scaler)

# Display image example
plt.figure(figsize=(16, 8))
im = plt.imshow(x_train[0], aspect='auto')
plt.colorbar(im)
plt.xlabel("Channel")
plt.ylabel("Sample / Position")
plt.title("EIS + DRT Input")
plt.xticks([0, 1, 2, 3],
    # [0,1], # EIS version only
    ["ZI", "ZII", "γ (gamma)","tau"]
    # ["ZI","ZII"] # EIS version only
)
plt.tight_layout()
plt.show()

# Output scaling
y_scaler = MinMaxScaler(feature_range=(0,1))

Y_train = y_scaler.fit_transform(y_train.reshape(-1,1))
Y_test = y_scaler.transform(y_test.reshape(-1,1))

# Print dataset shapes (ensure correct shapes are being passed to the model)
print("X_train:", X_train.shape)
print("Y_train:", Y_train.shape)

print("X_test:", X_test.shape)
print("Y_test:", Y_test.shape)

# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------

# Possible weighted MSE loss function if output is a vector containing SoC and SoH
# This would require changing the construct_feature function to return a vector y containing SoC and SoH
def weighted_mse(y_true, y_pred):
    error = tf.square(y_true - y_pred)
    loss_1 = tf.reduce_mean(error[:, 0])
    loss_2 = tf.reduce_mean(error[:, 1])
    loss = 0.9 * loss_1 + 0.1 * loss_2
    return loss

# CNN model
cnn = keras.models.Sequential()
cnn.add(keras.layers.Conv2D(128,kernel_size=(7,2),activation='relu',input_shape=X_train.shape[1:]))
cnn.add(keras.layers.MaxPooling2D(pool_size=(2,1)))
cnn.add(keras.layers.Conv2D(256,kernel_size=(7,1),activation='relu'))
cnn.add(keras.layers.MaxPooling2D(pool_size=(2,1)))
cnn.add(keras.layers.Flatten())
# kernel_regularizer=keras.regularizers.l2(0.001)
# cnn.add(keras.layers.GlobalAveragePooling2D())
cnn.add(keras.layers.Dense(128,activation='relu'))
cnn.add(keras.layers.Dense(64,activation='relu',))
# Regression output
cnn.add(keras.layers.Dense(1,activation='linear'))

optimizer = keras.optimizers.AdamW(learning_rate=0.001)

cnn.compile(loss='mse',optimizer=optimizer,metrics=['mae'])

try:
    history = cnn.fit(X_train,Y_train,validation_data=(X_test,Y_test),epochs=1000,batch_size=1024)

    plt.plot(history.history['mae'])
    plt.plot(history.history['val_mae'])

    plt.title('Mean Absolute Error')
    plt.ylabel('MAE')
    plt.xlabel('Epoch')
    plt.legend(['Train','Test'])

    plt.show()

except KeyboardInterrupt:
    print("Training interrupted. Saving model...")
    cnn.save("interrupted_model.keras")

# Save model / load model
cnn.save("cnn_model.keras") # Save model
# cnn = keras.models.load_model("cnn_model.keras") # Load model
cnn.summary()

# ---------------------------------------------------------------------------
# Post Training Analysis
# ---------------------------------------------------------------------------
def compute_rmse(y_true,y_pred):
    n = len(y_true)
    rmse = np.sqrt((1/n)*(y_true - y_pred)**2)
    return rmse

def compute_rmspe(y_true,y_pred):
    n = len(y_true)
    rmspe = np.sqrt((1/n)*((y_true-y_pred)/(y_true))**2)
    return rmspe

def compute_mea(y_true,y_pred):
    n = len(y_true)
    mea = (1/n)*abs(y_true-y_pred)
    return mea

# Process of creating fully new dataset
def preprocess_for_cnn(images, scalers):
    # Shape (N, 3, 510)
    images = np.stack(images)
    # Shape (N, 510, 3)
    images = np.transpose(images, (0, 2, 1))
    # Shape (N, 510, 3, 1)
    images = images[..., np.newaxis]
    # Apply training scalers
    for channel in range(2):
        images[:, :, channel, 0] = scalers[channel].transform(
            images[:, :, channel, 0]
        )
    return images

for c, drt, s in zip(data_cell10, drt_cell10, soc_cell10):
    # Construct dataset
    images_c, _ = construct_feature(directory_cell10,c,drt,0,0)

    # Preprocessing
    X_final = preprocess_for_cnn(images_c,scalers)

    # Prediction
    pred_scaled = cnn.predict(X_final,verbose=0)

    # Transform back to SOC
    pred = y_scaler.inverse_transform(pred_scaled)

    pred_soc = pred[:, 0]

    print(
        f"True SOC: {s:.2f}% | "
        f"Predicted SOC: {pred_soc.mean():.2f}% "
        f"+/- {pred_soc.std():.2f}% | "
        f"Difference: {((s-pred_soc.mean())/s)*100:.2f}% "
        f"True SOH: {soh_cell10:.2f}% | "
    )

# Prediction on out-of-sample data but for batteries that have been seen
pred_scaled = cnn.predict(X_test,verbose=0)
soc_list = soc_test_cell12 + soc_test_cell1 + soc_test_cell01

# Transform back to SOC 
pred = y_scaler.inverse_transform(pred_scaled)

pred_soc = pred[:,0]

diff = ((np.array(soc_list) - pred_soc.mean()) / np.array(soc_list)) * 100

for soc,pred,diff in zip(soc_list,pred_soc,diff):
    print(f"Actual SoC: {soc}, Predicted SoC: {pred.mean()}, Difference (%): {diff}")
