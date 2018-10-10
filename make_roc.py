import pandas as pd
import matplotlib.pyplot as plt

sfp_roc_data = pd.read_csv('sfp_specific_roc_data.csv')
plt.figure()
plt.plot(sfp_roc_data['x'], sfp_roc_data['y'])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Sfp')
plt.savefig('sfp_roc.pdf')


acps_roc_data = pd.read_csv('acps_specific_roc_data.csv')
plt.figure()
plt.plot(acps_roc_data['x'], acps_roc_data['y'])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Sfp')
plt.savefig('acps_roc.pdf')
