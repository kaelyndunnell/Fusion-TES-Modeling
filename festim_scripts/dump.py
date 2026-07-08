# train model 
from autoemulate.emulators import GaussianProcessRBF

def make_gp(x_train, y_train, lr=5e-2):
    return GaussianProcessRBF(
        x_train, 
        y_train, 
        lr=lr, 
        # standardize_x=True,
        standardize_y=False,
    )

# train initial emulator 
from sklearn.preprocessing import StandardScaler

# Standardize inputs
scaler = StandardScaler()
X_standardized = scaler.fit_transform(X)

# Train emulator on standardized inputs
print("training initial emulator...")
emulator = make_gp(X_standardized, Y, 0.1)
emulator.fit(X_standardized, Y)

# Save the scaler for later use
# dump(scaler, os.path.join(path, "input_scaler.joblib"))
# print("training initial emulator...")
# emulator = make_gp(X,Y, 0.1)
# emulator.fit(X,Y)
print("emulator trained!")
print("Lengthscale:", emulator.covar_module.base_kernel.lengthscale)

print("Input ranges (raw):")
for k, v in simulator.parameters_range.items():
    print(f"  {k}: {v}, span = {v[1] - v[0]}")

print("\nInput stats (after loading):")
print(f"  X shape: {X.shape}")
print(f"  X mean: {X.mean(axis=0)}")
print(f"  X std: {X.std(axis=0)}")

training initial emulator...
emulator trained!
Lengthscale: tensor([[[1.6611, 3.6127, 3.5336, 3.6474]],

        [[0.7413, 3.4765, 3.3059, 3.3481]]], grad_fn=<SoftplusBackward0>)
Input ranges (raw):
  bend_radius: (0.02, 0.2), span = 0.18000000000000002
  length: (0.5, 1.1), span = 0.6000000000000001
  v_in: (0.3114675730876479, 0.97), span = 0.658532426912352
  c_in: (np.float64(-7.779747962914968), np.float64(1.2202520370850323)), span = 9.0

Input stats (after loading):
  X shape: (74, 4)
  X mean: [-3.65602662  0.66445946  0.10608108  0.80189189]
  X std: [2.61824602 0.18965931 0.05153895 0.17620378]