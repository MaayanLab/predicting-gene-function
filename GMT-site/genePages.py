from flask import Flask, request, url_for, redirect, render_template
from wtforms import Form, StringField, validators
from string import Template
app = Flask(__name__)

HTML_TEMPLATE = Template("""
	<div class="topnav" id="myTopnav">
		<a href="{{ url_for('testhome') }}">Home</a>
	</div> 
	<h2>${gene_name}</h2>
	""")




class MyForm(Form):
    # username = StringField('Username', [validators.Length(min=4, max=25)])
    # email = StringField('Email Address', [validators.Length(min=6, max=35)])
    # password = PasswordField('New Password', [
    #     validators.DataRequired(),
    #     validators.EqualTo('confirm', message='Passwords must match')
    # ])
    # confirm = PasswordField('Repeat Password')
    # accept_tos = BooleanField('I accept the TOS', [validators.DataRequired()])

    gene = StringField('Gene Name', [validators.Length(min=1, max=25)])




@app.route('/register', methods=['GET', 'POST'])
def register():
    form = MyForm(request.form)
    if request.method == 'POST' and form.validate():
        return redirect(url_for('some_gene_page, some_gene = form.gene.data'))
    return render_template('register.html', form=form)





@app.route('/', methods = ['GET', 'POST'])
def testhome():
	form = MyForm(request.form)
	if request.method == 'POST' and form.validate():
		return redirect(url_for('some_gene_page', some_gene=form.gene.data))

	return render_template('test-home.html', form = form)


@app.route('/<some_gene>', methods = ['GET', 'POST'])
def some_gene_page(some_gene):
	if request.method == 'POST':
		return redirect(url_for('testhome'))
    	# do stuff when the form is submitted

    	# redirect to end the POST handling
    	# the redirect can be to the same route or somewhere else

	# show the form, it wasn't submitted

	with open('static/predicted_z_' + some_gene.upper() + '.txt') as gene_textFile:
		gene_list = gene_textFile.readlines()

	gene_list = [x.strip() for x in gene_list]

	# gene_list_split = []

	# gene_list_split = [line.split('\t') for line in gene_list]

	return render_template('some-gene.html', some_gene = some_gene, gene_list = gene_list)


@app.route('/hello/')
@app.route('/hello/<name>')
def hello(name=None):
	return render_template('hello.html', name = name)


@app.route('/history')
def history():

	with open('events.txt') as event_log:
		event_lines = event_log.readlines()

	door_history = []
	events = len(event_lines)

	for x in range(events):
		door_history.append(event_lines[x].split(" "))

	return render_template('door-log.html', events=events, door_history=door_history)


if __name__ == '__main__':
    	app.run(debug=True)