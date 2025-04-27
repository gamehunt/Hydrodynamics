from tabulate import tabulate

if __name__ == '__main__':
    indexes = []
    for i in range(3):
        for j in range(7):
            if i == 0 and j > 5:
                continue
            indexes.append(f'{j + 1}{i + 1}')
    matrix = [[''] * 22]
    matrix[0][20] = 'V'
    matrix[0][21] = 'F'
    with open('output.txt', 'w') as f:
        for i in range(20):
            row     = ['0'] * 22
            row[20] = 'V' + indexes[i]
            row[21] = '-f' + indexes[i]

            if i == 5 or i == 12:
                # плохие узлы
                if i == 5:
                    row[i]  = 'B'
                    row[i - 1] = '2/h^2(1 + a)'
                    row[i + 6] = '1/k^2'
                    row[21] += '-mu3(x)/k^2 - 2/(h^2a(1 - a))'
                else:
                    row[i]  = 'C'
                    row[i - 1] = '1/h^2'
                    row[i + 7] = '2/k^2(1 + b)'
                    row[21] += '-mu2(x)/h^2 - 2/(k^2b(1 + b))'
            else:
                row[i]  = 'A'
                # Левый сосед
                if i == 0 or i == 6 or i == 13:
                    row[21] += '-mu1(y)/h^2'
                else:
                    row[i - 1] = '1/h^2'

                # Правый сосед
                if i == 12 or i == 19:
                    row[21] += '-mu2(y)/h^2'
                else:
                    row[i + 1] = '1/h^2'

                # Нижний сосед
                if i < 5:
                    row[21] += '-mu3(x)/k^2'
                elif i == 12:
                    pass
                elif i > 5 and i < 13:
                    row[i - 6] = '1/k^2'
                else:
                    row[i - 7] = '1/k^2'

                # Верхний сосед
                if i > 12:
                    row[21] += '-mu4(x)/k^2'
                elif i < 6:
                    row[i + 6] = '1/k^2'
                else:
                    row[i + 7] = '1/k^2'

            matrix.append(row)
        f.write(tabulate(matrix, headers="firstrow"))
