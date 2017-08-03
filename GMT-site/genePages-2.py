# from flask import Flask, request, url_for, redirect, render_template

# app = Flask(__name__)

# @app.route('/')
# def index():
#     return render_template('tempIndex.html')

# @app.route('/cool_form', methods=['GET', 'POST'])
# def cool_form():
#     if request.method == 'POST':
#         # do stuff when the form is submitted

#         # redirect to end the POST handling
#         # the redirect can be to the same route or somewhere else
#         return redirect(url_for('index'))

#     # show the form, it wasn't submitted
#     return render_template('cool_form.html')

# if __name__ == '__main__':
# 	app.run(debug=True)


from flask import Flask, request, url_for, redirect, render_template
from string import Template
app = Flask(__name__)

HTML_TEMPLATE = Template("""
	<h1>${gene_name}</h1>
	""")

@app.route('/')
def home():
	return render_template('home(2).html')

# @app.route('/prkag1', methods=['GET', 'POST'])
# def prkag1():
# 	if request.method == 'POST':
#     # do stuff when the form is submitted

#     # redirect to end the POST handling
#     # the redirect can be to the same route or somewhere else
#     return redirect(url_for('home'))

# # show the form, it wasn't submitted
# return render_template('prkag1.html')

@app.route('/<some_gene>', methods = ['GET', 'POST'])
def some_gene_page(some_gene):
	if request.method == 'POST':
		return redirect(url_for('home'))
    	# do stuff when the form is submitted

    	# redirect to end the POST handling
    	# the redirect can be to the same route or somewhere else

	# show the form, it wasn't submitted
	return(HTML_TEMPLATE.substitute(gene_name=some_gene))

if __name__ == '__main__':
    	app.run(debug=True)