
import sys
import os
import pandas as pd
import numpy as np
import re

def drop_lines(text, linenos):
	lines = text.splitlines()
	linenos.sort(reverse=True)

	for lineno in linenos:
		del lines[lineno]
	return '\n'.join(lines)

def replace_line(text, lineno, newline):
	lines = text.splitlines()
	lines[lineno] = newline
	return '\n'.join(lines)

def nlines(text):
	lines = text.splitlines()
	return len(lines)

def apply_float_formatting(df):
	for col in df.columns:
		if col != 'decimals':
			new_col = df.apply(lambda x: float2string(x, col), axis=1)
			df[col] = new_col

	del df['decimals']

	return df

def float2string(data, colname):
	precision = int(data['decimals'])
	try:
		if np.isnan(data[colname]):
			return ''
		else:
			return f'{data[colname]:.{int(data["decimals"])}f}'
	except:
		return data[colname]

def header_panel(filepath):
	df = pd.read_excel(filepath, index_col=0, header=0)
	df = apply_float_formatting(df)
	colfmt = 'l'
	for ic in range(len(df.columns)):
		colfmt += 'c'
	tex = df.to_latex(float_format="%.1f", escape=False, na_rep='', column_format=colfmt)

	n = nlines(tex)
	lines_to_drop = [3, n-2, n-1]
	tex = drop_lines(tex, lines_to_drop)
	tex = re.sub(r'__v\d+__', '', tex)

	return tex

def other_panel(dirpath, table, panel, panelname):
	filename = f'table{table}_panel{panel}.xlsx'
	filepath = os.path.join(dirpath, filename)

	df = pd.read_excel(filepath, index_col=0, header=0)
	df = apply_float_formatting(df)
	colfmt = 'l'
	for ic in range(len(df.columns)):
		colfmt += 'c'
	tex = df.to_latex(float_format="%.1f", escape=False, na_rep='', column_format=colfmt)

	cols = len(df.columns)
	newline = ''
	for i in range(cols):
		newline += ' & '

	newline += ' \\\\'

	n = nlines(tex)
	tex = replace_line(tex, 0, newline)
	tex = drop_lines(tex, [3, n-2, n-1])

	headername = f'Panel {panel}: {panelname}'
	newline = f'\\multicolumn{{{cols+1}}}{{c}}{{\\textbf{{{headername}}}}}\\\\'
	tex = replace_line(tex, 2, newline)
	tex = re.sub(r'__v\d+__', '', tex)

	return tex

def save_tex_table(dirpath, tableno):
	panels = [header_panel(os.path.join(dirpath, f'table{tableno}_header.xlsx'))]

	if tableno == 1:
		panel_args = [
			['A', 'Income Statistics'],
			['B', 'Wealth Statistics'],
			['C', 'MPC Size Effects'],
			['D', 'MPC Sign Effects'],
		]
	elif tableno == 2:
		panel_args = [
			['A', 'Decomposition of Mean MPC, around 0'],
			['B', 'Decomposition of Mean MPC, around 0.01'],
			['C', 'Decomposition of Mean MPC, around 0.05'],
			['D', 'Decomposition of Mean MPC - MPC$_{RA}$'],
		]
	else:
		panel_args = [
			['A', 'Quarterly MPC Decomp w.r.t. Baseline'],
			['A2', 'Quarterly MPC Decomp as \\% of $E[m_1]-E[m_0]$'],
			['B', 'Wealth Statistics'],
			['C', 'MPC Size Effects'],
			['D', 'MPC Sign Effects'],
		]

	for ip in range(len(panel_args)):
		panels.append(other_panel(dirpath, tableno, *panel_args[ip]))

	tex = '\n'.join(panels)
	tex += '\n\\bottomrule\n\\end{tabular}'
	
	texfilepath = os.path.join(dirpath, f'table{tableno}.tex')
	with open(texfilepath, 'w') as fobj:
		fobj.write(tex)

if __name__ == '__main__':
	dirpath = sys.argv[1]

	for it in range(1, 10):
		save_tex_table(dirpath, it)