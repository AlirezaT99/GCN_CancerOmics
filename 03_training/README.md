## Train Mogonet

Having created the input matrices, we should feed the data to Mogonet architecture. Following these steps:

1. Clone the repository: (preferably in Colab or a server with GPU available)

```bash
git clone https://github.com/txWang/MOGONET.git
```

2. Configure the code. To do so, copy the output from previous stage to the `MOGONET` directory and add the directory name and info to `main_biomarker.py` and `main_mogonet.py`.

3. Run `main_mogonet.py` to train the network and store the model

4. Run `main_biomarker.py` to evaluate the feature importances and extract biomarkers from the model.

### Alternative approach

Or you can configurate the `runner.py` code and run it.