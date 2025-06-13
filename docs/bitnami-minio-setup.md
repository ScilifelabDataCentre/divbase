# Steps Taken to Create a MinIO Deployment Using Bitnami's Helm Chart

## Local Deployment Using k3d

```bash
# create a local k3d cluster
k3d cluster create mycluster

kustomize build kustomize/overlays/local --enable-helm | kubectl apply -f -

# access console - see values-local.yaml for username and password. 
kubectl port-forward service/minio 9001:9001 -n divbase-local

# make API available 
k port-forward service/minio 9000:9000 -n divbase-local
```


## Remote deployment on scilifelab-2-dev.sys.kth.se

Namespace `divbase-testground` already created via Rancher GUI.

### Step 1: Create secret for minio admin credentials - one time operation

```bash
# password removed from here but can be obtained from the rancher console
k --namespace divbase-testground \
    create secret generic minio-credentials \
    --dry-run=client \
    --from-literal=MINIO_ROOT_USER="divbaseadmin" \
    --from-literal=MINIO_ROOT_PASSWORD="" \ 
    -o json \
    | kubeseal -o yaml > charts/bitnami-minio/secrets/secret-minio-seal.yaml

k apply -f charts/bitnami-minio/secrets/secret-minio-seal.yaml 
```

### Step 2: kustomize build and apply

```bash
kustomize build kustomize/overlays/scilifelab-2-dev --enable-helm | kubectl apply -f -
```