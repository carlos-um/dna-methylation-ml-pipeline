import pandas as pd
from sklearn.model_selection import train_test_split

df = pd.read_csv("../data/LBDfcFilteredMyNorm.csv", index_col=0)

# For testing
# df = df.sample(n=10000, random_state=123)

df = df.T

if "Unnamed: 0" in df.index:
    df = df.drop("Unnamed: 0", axis=0)

def extraer_categoria(nombre):
    if nombre.startswith("CTRL"):
        return "CTRL"
    elif nombre.startswith("PDD"):
        return "PDD"
    elif nombre.startswith("PD"):
        return "PD"
    elif nombre.startswith("DLB"):
        return "DLB"
    else:
        return None

df['Categoria'] = df.index.to_series().apply(extraer_categoria)

if df['Categoria'].isnull().any():
    print("Algunos nombres de fila no pudieron clasificarse en una categoría.")

X = df.drop("Categoria", axis=1)
y = df["Categoria"]

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.3, stratify=y, random_state=123
)

print("Índices del conjunto de entrenamiento:")
print(X_train.index.tolist())
print("\nÍndices del conjunto de prueba:")
print(X_test.index.tolist())

X_train = X_train.T

X_train.to_csv("../data/LBDfcFilteredMyNorm_train.csv")

X_test = X_test.T

X_test.to_csv("../data/LBDfcFilteredMyNorm_test.csv")

# Sheet

df_sheet = pd.read_csv("../data/LBDfcSamplesheet.csv")

train_indices = X_train.columns.tolist()
test_indices = X_test.columns.tolist()

df_sheet_train = df_sheet[df_sheet["Sample_Name"].isin(train_indices)]
df_sheet_train = df_sheet_train.set_index("Sample_Name").reindex(X_train.columns).reset_index()

df_sheet_test = df_sheet[df_sheet["Sample_Name"].isin(test_indices)]
df_sheet_test = df_sheet_test.set_index("Sample_Name").reindex(X_test.columns).reset_index()

df_sheet_train.to_csv("../data/LBDfcSamplesheet_train.csv", index=False)
df_sheet_test.to_csv("../data/LBDfcSamplesheet_test.csv", index=False)
