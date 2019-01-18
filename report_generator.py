import copy
import datetime
import random
import itertools
import numpy as np
from pylatex import Document, Section, Subsection, LongTable, Matrix, Math, Itemize, Figure, NoEscape
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
class Model:
    def __init__(self, model_name, model_info, training_metrics, testing_metrics):
        self.model_name = model_name
        self.model_info = model_info
        self.training_metrics = training_metrics
        self.testing_metrics = testing_metrics

    def get_model_name_and_info(self):
        return self.model_name, self.model_info

    def get_model_training_testing_metrics(self):
        return self.training_metrics, self.testing_metrics

def fill_initial_info(doc, dataset_name, num_indivs, n_snps,training_shape, test_shape, k_folds):

    with doc.create(Itemize()) as itemize:
        itemize.add_item('Dataset Name: {}'.format(dataset_name))
        itemize.add_item('# Individuals: {}'.format(num_indivs))
        itemize.add_item('# SNPs: {}'.format(n_snps))
        tot = (training_shape[0] + test_shape[0]) / 100
        itemize.add_item('Train-Test split: {}-{}'.format(training_shape[0] / tot, test_shape[0] / tot))
        itemize.add_item('{}-fold Cross Validation'.format(k_folds))

def fill_info(doc, method_info):
    with doc.create(Itemize()) as itemize:
        for info in method_info:
            itemize.add_item(info)

def fill_metrics_table(doc, metrics, col_string, models, index):
    copy_metrics = copy.deepcopy(metrics)
    copy_metrics.insert(0, 'Method')
    with doc.create(LongTable(col_string)) as data_table:
        data_table.add_hline()
        data_table.add_row(copy_metrics)
        data_table.add_hline()
        data_table.end_table_header()
        data_table.add_hline()
        for model in models:
            to_show = model.get_model_training_testing_metrics()[index]
            to_show.insert(0, model.model_name)
            data_table.add_row(to_show)


	
def fill_confusion_matrix(doc, cm):
	matrix = Matrix(cm, mtype='b')
	math = Math(data=['M=', matrix])
	doc.append(math)
	accuracy = np.trace(cm) / float(np.sum(cm))
	misclass = 1 - accuracy
	target_names = ['Controls','Cases']
	cmap = plt.get_cmap('Blues')
	fig = plt.figure(figsize=(10, 7))
	plt.imshow(cm, interpolation='nearest', cmap=cmap)
	plt.title('Confusion Matrix',fontsize=18)
	plt.colorbar()
	normalize = True
	if target_names is not None:
		tick_marks = np.arange(len(target_names))
		plt.xticks(tick_marks, target_names, rotation=45)
		plt.yticks(tick_marks, target_names)

	if normalize:
		cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]


	thresh = cm.max() / 1.5 if normalize else cm.max() / 2
	for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
		if normalize:
			plt.text(j, i, "{0:0.4f}".format(cm[i, j]),
					 fontsize=18,
					 horizontalalignment="center",
					 color="white" if cm[i, j] > thresh else "black")
		else:
			plt.text(j, i, "{:,}".format(cm[i, j]),
					fontsize=18,
					horizontalalignment="center",
					color="white" if cm[i, j] > thresh else "black")


	plt.tight_layout()
	plt.ylabel('True label',fontsize=18)
	plt.xlabel('Predicted label\naccuracy={:0.4f}; misclass={:0.4f}'.format(accuracy, misclass),fontsize=18)
	#plt.show()
	with doc.create(Figure(position='htbp')) as plot:
		plot.add_plot()
		plot.add_caption('Confusion Matrix')

def fill_document(doc, models, metrics, confusion_matrix, dataset_name, num_indivs, n_snps,training_shape, test_shape, k_folds):

    """
    The main method that does the document creation in LaTeX
    :param metrics: The list of metrics' names that were stored in the models. Must be in same order.
    :param models: List of Model objects containing relevant info per model.
    :param k_folds: number of folds in k-fold cv
    :param training_shape: numpy shape as (rows, cols)
    :param test_shape: numpy shape as (rows, cols)
    :param n_snps: number of SNPs
    :param num_indivs: number of individuals
    :param dataset_name: name of dataset
    :param doc: The doc file. This is the python representation of our file
    :param confusion_matrix: a NUMPY array containing the confusion matrix.
    """
    with doc.create(Section('Data Info:')):
        fill_initial_info(doc, dataset_name, num_indivs, n_snps,
                          training_shape, test_shape, k_folds)
    with doc.create(Section('Models')):
        for model in models:
            model_name, model_info = model.get_model_name_and_info()
            with doc.create(Subsection('{}'.format(model_name))):
                fill_info(doc, model_info)
    num_cols = len(models[0].get_model_training_testing_metrics()[0])
    col_string = " ".join(["l"] * (num_cols + 1))
    with doc.create(Section('Metrics: Training')):
        fill_metrics_table(doc, metrics, col_string, models, 0)
    with doc.create(Section('Metrics: Testing')):
        fill_metrics_table(doc, metrics, col_string, models, 1)
    with doc.create(Section('Confusion Matrix')):
        fill_confusion_matrix(doc, confusion_matrix)
    doc.generate_pdf(clean_tex=True)
    doc.generate_tex()

def generate_pdf_with_name(prefix,rn):
    # time_string = datetime.datetime.now()
    # rn = random.randint(1,200000)
    file_name = "{}_{}".format(prefix, rn)
    doc = Document(file_name)
    return doc

if __name__ == '__main__':
    # Basic document
    doc = Document('basic66')
    M = np.matrix([[2, 3, 4],
                   [0, 0, 1],
                   [0, 0, 2]])
    model1 = Model('Model1', ['SomeHyperParam=3'], [0.78, 0.9], [0.75, 0.8])
    model2 = Model('Model2', ['SomeHyperParam=3', 'OtherHyperParam=0.9'], [0.98, 0.5], [0.45, 0.4])
    model3 = Model('Model3', [], [0.88, 0.6], [0.85, 0.4])
    model4 = Model('Model4', ['SomeHyperParam=3'], [0.6, 0.7], [0.62, 0.67])

    # fill_document(doc, models, metrics, confusion_matrix, dataset_name, num_indivs, n_snps,training_shape, test_shape, k_folds):
    fill_document(doc, [model1, model2, model3, model4], ['acc', 'f1'], M, 'Some_Dataset', 1000, 12000, (700, 12000),
                  (300, 12000), 10)
    doc.generate_pdf(clean_tex=True)
    doc.generate_tex()
    tex = doc.dumps()  # The document as string in LaTeX syntax
