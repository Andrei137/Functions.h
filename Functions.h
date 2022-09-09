#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cstring>
#include <ctime>
#include <cstring>
#include <windows.h>
#include <conio.h>


class generalPurpose_ALGO
{
    public :
        // O(firstLength + secondLength)
        template <typename arr> inline
        void merge(int firstLength, arr firstArray[], int secondLength, arr secondArray[], int &mergedLength, arr mergedArray[])
        {
            int i = 0, j = 0;
            mergedLength = 0;
            while (i < firstLength && j < secondLength)
                if (firstArray[i] < secondArray[j])
                    mergedArray[mergedLength++] = firstArray[i++];
                else
                    mergedArray[mergedLength++] = secondArray[j++];
            while (i < firstLength)
                mergedArray[mergedLength++] = firstArray[i++];
            while (j < secondLength)
                mergedArray[mergedLength++] = secondArray[j++];
        }

        // O(log(number))
        int numberOfDigits(int number)
        {
            int nrDigits = 0;
            do
                ++nrDigits;
            while (number /= 10);
            return nrDigits;
        }

        // O(sqrt(number))
        int numberOfPrimeFactors(int number)
        {
            int div = 2, nrPrimeDiv = 0;
            while (number != 1)
            {
                if (number % div == 0)
                {
                    while (number % div == 0)
                        number /= div;
                    ++nrPrimeDiv;
                }
                ++div;
                if (number > 1 && div * div > number)
                    div = number;
            }
            return nrPrimeDiv;
        }

        // O(log(number))
        int reversedNumber(int number)
        {
            int numberReveresed = 0;
            do
                numberReveresed = numberReveresed * 10 + number % 10;
            while (number /= 10);
            return numberReveresed;
        }

};


class divsAndMx_ALGO
{
    public :
        // O(log(min(firstNumber, secondNumber))
        int greatestCommonDivisor(int firstNumber, int secondNumber)
        {
            int remainder;
            do
            {
                remainder = firstNumber % secondNumber;
                firstNumber = secondNumber;
                secondNumber = remainder;
            } while (secondNumber);
            return firstNumber;
        }

        // O(log(min(firstNumber, secondNumber))
        int leastCommonMultiple(int firstNumber, int secondNumber)
        {
            return firstNumber * secondNumber / greatestCommonDivisor(firstNumber, secondNumber);
        }

        // O(sqrt(number))
        int numberOfDivisors(int number)
        {
            int div, nrDiv = 2;
            for (div = 2; div * div < number; ++div)
                if (number % div == 0)
                    nrDiv += 2;
            if (div * div == number)
                ++nrDiv;
            return nrDiv;
        }
};


class primality_ALGO
{
    public :
        // O(sqrt(number))
        bool isPrime(int number)
        {
            if (number == 2)
                return true;
            if (number < 2 || number % 2 == 0)
                return false;
            int div = 3;
            while (number % div && div * div < number)
                div += 2;
            return div * div > number;
        }

        // O(length * log(log(length)))
        template <typename arr> inline
        void sieveOfEratosthenes(int length, arr &isPrime)
        {
            int i, j;
            isPrime[0]  = isPrime[1] = false;
            for (i = 2; i < length; ++i)
                isPrime[i] = true;
            for (i = 2; i * i <= length; ++i)
                if (isPrime[i])
                    for (j = i; j * i <= length; ++j)
                        isPrime[i * j] = false;
        }
};


class power_ALGO
{
    public :
        // O(log2(power))
        int fastExponentiationIterative(int number, int power)
        {
            int powerOfNumber = 1;
            do
            {
                if (power % 2)
                    powerOfNumber *= number;
                number *= number;
            } while (power /= 2);
            return powerOfNumber;
        }

        // O(log2(power))
        int fastExponentiationRecursive(int number, int power)
        {
            if (!power)
                return 1;
            if (power % 2)
                return number * fastExponentiationRecursive(number, power - 1);
            int powerOfNumber = fastExponentiationRecursive(number, power / 2);
            return powerOfNumber * powerOfNumber;
        }

        // O(power)
        int power(int number, int power)
        {
            int powerOfNumber = 1;
            while (power--)
                powerOfNumber *= number;
            return powerOfNumber;
        }
};


// o(log(number))
class conversion_ALGO
{
    private :
        #define elif else if

        // Used for convertRomanToInt()
        bool charEquals(char c, char x, char y)
        {
            return (c == x || c == y);
        }

        // Used for convertRomanToInt()
        int charToInt(char c)
        {
            if (c == 'I')
                return 1;
            elif (c == 'V')
                 return 5;
            elif (c == 'X')
                return 10;
            elif (c == 'L')
                return 50;
            elif (c == 'C')
                return 100;
            elif (c == 'D')
                return 500;
            return 1000;
        }
    public :
        int convert10Tob(int number, int base)
        {
            int convertedNumber = 0, power = 1;
            do
            {
                convertedNumber += power * (number % base);
                power *= 10;
            }
            while (number /= base);
            return convertedNumber;
        }

        int convertbTo10(int number, int base)
        {
            int convertedNumber = 0, power = 1;
            do
            {
                convertedNumber += number % 10 * power;
                power *= base;
            }
            while (number /= 10);
            return convertedNumber;
        }

        int convertRomanToInt(std::string s)
        {
            int n = 0, i, size = s.size() - 1;
            for (i = 0; i < size; ++i)
                if (s[i] == 'I' && charEquals(s[i + 1], 'V', 'X'))
                    --n;
                elif (s[i] == 'X' && charEquals(s[i + 1], 'L', 'C'))
                    n -= 10;
                elif (s[i] == 'C' && charEquals(s[i + 1], 'D', 'M'))
                    n -= 100;
                else
                    n += charToInt(s[i]);
            return n + charToInt(s[size]);
        }
};


// O(arrayLength)
class blockSwap_ALGO
{
    public :
        template <typename arr> inline
        void blockSwapIterative(int arrayLength, arr &array, int section)
        {
            int i, j, k, temp;
            if (!section || section == arrayLength)
                return;
            if(section > arrayLength)
              section %= arrayLength;
            i = section;
            j = arrayLength - section;
            while (i != j)
            {
                if (i < j)
                {
                    for (k = 0; k < i; ++k)
                    {
                        temp = array[section - i + k];
                        array[section - i + k] = array[section + j - i + k];
                        array[section + j - i + k] = temp;
                    }
                    j -= i;
                }
                else
                {
                    for (k = 0; k < j; ++k)
                    {
                        temp = array[arrayLength - i + k];
                        array[arrayLength - i + k] = array[arrayLength + k];
                        array[arrayLength + k] = temp;
                    }
                    i -= j;
                }
            }
            for (k = 0; k < i; ++k)
            {
                temp = array[arrayLength - i + k];
                array[arrayLength - i + k] = array[arrayLength + k];
                array[arrayLength + k] = temp;
            }
        }

        template <typename arr> inline
        void blockSwapRecursive(int arrayLength, arr &array, int section)
        {
            int i, temp;
            if(!section || section == arrayLength)
                return;
            if(section > arrayLength)
              section %= arrayLength;
            if(arrayLength - section == section)
            {
                for(i = 0; i < section; ++i)
                {
                    temp = array[i];
                    array[i] = array[arrayLength - section + i];
                    array[arrayLength - section + i] = temp;
                }
                return;
            }
            if(section < arrayLength - section)
            {
                for(i = 0; i < section; i++)
                {
                    temp = array[i];
                    array[i] = array[arrayLength - section + i];
                    array[arrayLength - section + i] = temp;
                }
                blockSwapRecursive(arrayLength - section, array, section);
            }
            else
            {
                for(i = 0; i < arrayLength - section; i++)
                {
                    temp = array[i];
                    array[i] = array[section + i];
                    array[section + i] = temp;
                }
                blockSwapRecursive(section, array + arrayLength - section, 2 * section - arrayLength);
            }
        }
};


class searching_ALGO
{
    public :
        // O(log(arrayLength))
    template <typename arr> inline
        int binarySearchIterative(int arrayLength, arr &array, int wantedNumber)
        {
            int left = 0, right = arrayLength, middle;
            while (left < right)
            {
                middle = (left + right) / 2;
                if (array[middle] == wantedNumber)
                    return middle;
                else if (array[middle] < wantedNumber)
                    left = middle + 1;
                else
                    right = middle - 1;
            }
            return -1;
        }

        // O(log(arrayLength))
        template <typename arr> inline
        int binarySearchRecursive(arr &array, int left, int right, int wantedNumber)
        {
            if (left <= right)
            {
                int middle = left + (right - left) / 2;
                if (array[middle] == wantedNumber)
                    return middle;
                if (array[middle] > wantedNumber)
                    return binarySearchRecursive(array, left, middle - 1, wantedNumber);
                return binarySearchRecursive(array, middle + 1, right, wantedNumber);
            }
            return -1;
        }

        // O(log(arrayLength))
        template <typename arr> inline
        int exponentialSearch(int arrayLength, arr &array, int wantedNumber)
        {
            if (array[0] == wantedNumber)
                return 0;
            int i = 1, minimum;
            while (i < arrayLength && array[i] <= wantedNumber)
                i *= 2;
            minimum = (i < arrayLength - 1) ? i : (arrayLength - 1);
            return binarySearchRecursive(array, i / 2, minimum, wantedNumber);
        }

        // O(log(log(arrayLength)))
        template <typename arr> inline
        int interpolationSearch(int arrayLength, arr &array, int wantedNumber)
        {
            int left = 0, right = arrayLength - 1, position;
            while (left <= right && wantedNumber >= array[left] && wantedNumber <= array[right])
            {
                if (left == right)
                {
                    if (array[left] == wantedNumber)
                        return left;
                    return -1;
                }
                position = left + (((double)(right - left) / (array[right] - array[left])) * (wantedNumber - array[left]));
                if (array[position] == wantedNumber)
                    return position;
                if (array[position] < wantedNumber)
                    left = position + 1;
                else
                    right = position - 1;
            }
            return -1;
        }

        // O(sqrt(arrayLength))
        template <typename arr> inline
        int jumpSearch(int arrayLength, arr &array, int wantedNumber)
        {
            int step, position = 0, minimum;
            //step = sqrt(arrayLength);
            minimum = (step < arrayLength) ? step : arrayLength;
            while (array[minimum - 1] < wantedNumber)
            {
                position = step;
                //step += sqrt(arrayLength);
                minimum = (step < arrayLength) ? step : arrayLength;
                if (position >= arrayLength)
                    return -1;
            }
            while (array[position] < wantedNumber)
            {
                ++position;
                minimum = (step < arrayLength) ? step : arrayLength;
                if (position == minimum)
                    return -1;
            }
            if (array[position] == wantedNumber)
                return position;
            return -1;
        }

        // O(arrayLength)
        template <typename arr> inline
        int linearSearch(int arrayLength, arr &array, int wantedNumber)
        {
            int i;
            for (i = 0; i < arrayLength; ++i)
                if (array[i] == wantedNumber)
                    return i;
            return -1;
        }

        // O(log(arrayLength))
        template <typename arr> inline
        int ternarySearchIterative(arr &array, int left, int right, int wantedNumber)
        {
            while (left <= right)
            {
                int firstMiddle = left + (right - left) / 3, secondMiddle = right - (right - left) / 3;
                if (array[firstMiddle] == wantedNumber)
                    return firstMiddle;
                if (array[secondMiddle] == wantedNumber)
                    return secondMiddle;
                if (wantedNumber < array[firstMiddle])
                    right = firstMiddle - 1;
                else if (wantedNumber > array[secondMiddle])
                    left = secondMiddle + 1;
                else
                {
                    left = firstMiddle + 1;
                    right = secondMiddle - 1;
                }
            }
            return -1;
        }

        // O(log(arrayLength))
        template <typename arr> inline
        int ternarySearchRecursive(arr &array, int left, int right, int wantedNumber)
        {
            if (left <= right)
            {
                int firstMiddle = left + (right - left) / 3, secondMiddle = right - (right - left) / 3;
                if (array[firstMiddle] == wantedNumber)
                    return firstMiddle;
                if (array[secondMiddle] == wantedNumber)
                    return secondMiddle;
                if (wantedNumber < array[firstMiddle])
                    return ternarySearchRecursive(array, left, firstMiddle - 1, wantedNumber);
                else if (wantedNumber > array[secondMiddle])
                    return ternarySearchRecursive(array, secondMiddle + 1, right, wantedNumber);
                else
                    return ternarySearchRecursive(array, firstMiddle + 1, secondMiddle - 1, wantedNumber);
            }

            return -1;
        }
};


// O(arrayLength * log(arrayLength))
class fastSorting_ALGO
{
    private :
        // Used for inPlaceMergeSort()
        template <typename arr> inline
        void inPlaceMerge(arr &array, int firstStart, int middle, int end)
        {
            int secondStart = middle + 1;
            if (array[middle] <= array[secondStart])
                return;
            while (firstStart <= middle && secondStart <= end)
                if (array[firstStart] <= array[secondStart])
                    ++firstStart;
                else
                {
                    int value = array[secondStart], index = secondStart;
                    while (index != firstStart)
                    {
                        array[index] = array[index - 1];
                        --index;
                    }
                    array[firstStart] = value;
                    ++firstStart;
                    ++middle;
                    ++secondStart;
                }
        }

        // Used for heapSort()
        template <typename arr> inline
        void heapify(int arrayLength, arr &array, int i)
        {
            int largest = i, left = 2 * i + 1, right = 2 * i + 2, aux;
            if (left < arrayLength && array[left] > array[largest])
                largest = left;
            if (right < arrayLength && array[right] > array[largest])
                largest = right;
            if (largest != i)
            {
                aux = array[i];
                array[i] = array[largest];
                array[largest] = aux;
                heapify(arrayLength, array, largest);
            }
        }
    public :
        template <typename arr> inline
        void heapSort(int arrayLength, arr &array)
        {
            int i, aux;
            for (i = arrayLength / 2 - 1; i >= 0; i--)
                heapify(arrayLength, array, i);
            for (i = arrayLength - 1; i > 0; i--)
            {
                aux = array[0];
                array[0] = array[i];
                array[i] = aux;
                heapify(i, array, 0);
            }
        }

        template <typename arr> inline
        void inPlaceMergeSort(arr &array, int left, int right)
        {
            if (left < right)
            {
                int middle = left + (right - left) / 2;
                inPlaceMergeSort(array, left, middle);
                inPlaceMergeSort(array, middle + 1, right);
                inPlaceMerge(array, left, middle, right);
            }
        }

        template <typename arr> inline
        void mergeSort(arr &array, int left, int right, int temp[])
        {
            //temp[i] = 0 before input
            if (left < right)
            {
                int middle = (left + right) / 2;
                mergeSort(array, left, middle, temp);
                mergeSort(array, middle + 1, right, temp);
                int i = left, j = middle + 1, tempLength = 0;
                while (i <= middle && j <= right)
                    if (array[i] < array[j])
                        temp[++tempLength] = array[i++];
                    else
                        temp[++tempLength] = array[j++];
                while (i <= middle)
                    temp[++tempLength] = array[i++];
                while (j <= right)
                    temp[++tempLength] = array[j++];
                for (i = left, j = 1; i <= right; ++i, ++j)
                    array[i] = temp[j];
            }
        }

        template <typename arr> inline
        void quickSort(arr &array, int left, int right)
        {
            if (left < right)
            {
                int middle = (left + right) / 2, aux = array[left], i = left, j = right, d = 0;
                array[left] = array[middle];
                array[middle] = aux;
                while (i < j)
                {
                    if (array[i] > array[j])
                    {
                        aux = array[i];
                        array[i] = array[j];
                        array[j] = aux;
                        d = 1 - d;
                    }
                    i += d;
                    j -= 1 - d;
                }
                quickSort(array, left, i - 1);
                quickSort(array, i + 1, right);
            }
        }

        template <typename arr> inline
        void shellSort(int arrayLength, arr &array)
        {
            int gap, i, j, temp;
            for (gap = arrayLength / 2; gap > 0; gap /= 2)
            {
                for (i = gap; i < arrayLength; ++i)
                {
                    temp = array[i];
                    for (j = i; j >= gap && array[j - gap] > temp; j -= gap)
                        array[j] = array[j - gap];
                    array[j] = temp;
                }
            }
        }
};


// O(arrayLength^2)
class slowSorting_ALGO
{
    private :
        // Used for combSort()
        int combGap(int gap)
        {
            gap = (gap * 10) / 13;
            return gap < 1 ? 1 : gap;
        }
    public :
        template <typename arr> inline
        void bubbleSortIterative(int arrayLength, arr &array)
        {
            int i, aux;
            bool sorted;
            do
            {
                sorted = true;
                for (i = 0; i < arrayLength - 1; ++i)
                    if (array[i] > array[i + 1])
                    {
                        aux = array[i];
                        array[i] = array[i + 1];
                        array[i + 1] = aux;
                        sorted = false;
                    }
            }
            while (!sorted);
        }

        template <typename arr> inline
        void bubbleSortRecursive(int arrayLength, arr &array)
        {
            if (arrayLength == 1)
                return;
            int i, aux;
            for (i = 0; i < arrayLength - 1; ++i)
                if (array[i] > array[i + 1])
                {
                    aux = array[i];
                    array[i] = array[i + 1];
                    array[i + 1] = aux;
                }
            bubbleSortRecursive(arrayLength - 1, array);
        }

        template <typename arr> inline
        void cocktailSort(int arrayLength, arr &array)
        {
            bool swapped = true;
            int i, start = 0, end = arrayLength - 1, aux;
            while (swapped)
            {
                swapped = false;
                for (i = start; i < end; ++i)
                {
                    if (array[i] > array[i + 1])
                    {
                        aux = array[i];
                        array[i] = array[i + 1];
                        array[i + 1] = aux;
                        swapped = true;
                    }
                }
                if (!swapped)
                    break;
                swapped = false;
                --end;
                for (i = end - 1; i >= start; --i)
                {
                    if (array[i] > array[i + 1])
                    {
                        aux = array[i];
                        array[i] = array[i + 1];
                        array[i + 1] = aux;
                        swapped = true;
                    }
                }
                ++start;
            }
        }

        template <typename arr> inline
        void combSort(int arrayLength, arr &array)
        {
            int gap = arrayLength, i, aux;
            bool swapped = true;
            while (gap != 1 || swapped)
            {
                gap = combGap(gap);
                swapped = false;
                for (i = 0; i < arrayLength - gap; ++i)
                {
                    if (array[i] > array[i + gap])
                    {
                        aux = array[i];
                        array[i] = array[i + gap];
                        array[i + gap] = aux;
                        swapped = true;
                    }
                }
            }
        }

        template <typename arr> inline
        void cycleSort(int arrayLength, arr &array)
        {
            int writes = 0, cycleStart, i, item, position, aux;
            for (cycleStart = 0; cycleStart <= arrayLength - 2; ++cycleStart)
            {
                item = array[cycleStart];
                position = cycleStart;
                for (i = cycleStart + 1; i < arrayLength; ++i)
                    if (array[i] < item)
                        ++position;
                if (position == cycleStart)
                    continue;
                while (item == array[position])
                    ++position;
                if (position != cycleStart)
                {
                    aux = item;
                    item = array[position];
                    array[position] = aux;
                    ++writes;
                }
                while (position != cycleStart)
                {
                    position = cycleStart;
                    for (i = cycleStart + 1; i < arrayLength; ++i)
                        if (array[i] < item)
                            position += 1;
                    while (item == array[position])
                        position += 1;
                    if (item != array[position])
                    {
                        aux = item;
                        item = array[position];
                        array[position] = aux;
                        writes++;
                    }
                }
            }
        }

        template <typename arr> inline
        void gnomeSort(int arrayLength, arr &array)
        {
            int index = 0, aux;
            while (index < arrayLength)
            {
                if (!index)
                    ++index;
                if (array[index] >= array[index - 1])
                    ++index;
                else
                {
                    aux = array[index];
                    array[index] = array[index - 1];
                    array[index - 1] = aux;
                    --index;
                }
            }
        }

        template <typename arr> inline
        void insertionSortIterative(int arrayLength, arr &array)
        {
            int i, position, aux;
            for (i = 1; i < arrayLength; ++i)
            {
                position = i;
                while (position > 0 && array[position] < array[position - 1])
                {
                    aux = array[position];
                    array[position] = array[position - 1];
                    array[position - 1] = aux;
                    --position;
                }
            }
        }

        template <typename arr> inline
        void insertionSortRecursive(int arrayLength, arr &array)
        {
            if (arrayLength <= 1)
                return;
            insertionSortRecursive(arrayLength - 1, array);
            int last = array[arrayLength - 1], j = arrayLength - 2;
            while (j >= 0 && array[j] > last)
            {
                array[j + 1] = array[j];
                --j;
            }
            array[j + 1] = last;
        }

        template <typename arr> inline
        void oddEvenSort(int arrayLength, arr &array)
        {
            int i, aux;
            bool isSorted = false;
            while (!isSorted)
            {
                isSorted = true;
                for (i = 1; i <= arrayLength - 2; i += 2)
                    if (array[i] > array[i + 1])
                    {
                        aux = array[i];
                        array[i] = array[i + 1];
                        array[i + 1] = aux;
                        isSorted = false;
                    }
                for (i = 0; i <= arrayLength - 2; i += 2)
                    if (array[i] > array[i+1])
                    {
                        aux = array[i];
                        array[i] = array[i + 1];
                        array[i + 1] = aux;
                        isSorted = false;
                    }
            }
        }

        template <typename arr> inline
        void selectionSort(int arrayLength, arr &array)
        {
            int i, j, position, aux;
            for (i = 0; i < arrayLength - 1; ++i)
            {
                position = i;
                for (j = i + 1; j < arrayLength; ++j)
                    if (array[j] < array[position])
                        position = j;
                aux = array[i];
                array[i] = array[position];
                array[position] = aux;
            }
        }
};


class Ggraphs_ALGO
{
    private :
        #define INFINITY 1000000
        #define MAX 101

        struct Edge
        {
            int x, y, c;
        };

        // Used for dijkstraWithPathShow()
        template <typename arr> inline
        void dijkstraFindPath(int currentVertex, int startingVertex, arr &previousVertex, int &length, arr &path)
        {
            length = 0;
            if(currentVertex != startingVertex)
                dijkstraFindPath(previousVertex[currentVertex], startingVertex, previousVertex, length, path);
            path[length++] = currentVertex;
        }

        // Used for primWithParentsShow()
        template <typename arr> inline
        void primWithParents(int nrVertices, int nrEdges, int parentVertex, int &cost, Edge edge[], arr &parent, arr &selected)
        {
            int i, j;
            for (i = 1; i <= nrVertices; ++i)
                parent[i] = selected[i] = 0;
            cost = 0;
            selected[parentVertex] = true;
            quickSortGraphs(edge, 0, nrEdges - 1);
            for (i = 1; i < nrVertices; ++i)
            {
                j = 0;
                while (selected[edge[j].x] == selected[edge[j].y] && j < nrEdges)
                    ++j;
                cost += edge[j].c;
                if (!selected[edge[j].x])
                {
                    selected[edge[j].x] = true;
                    parent[edge[j].x] = edge[j].y;
                }
                else
                {
                    selected[edge[j].y] = true;
                    parent[edge[j].y] = edge[j].x;
                }
            }
        }

        // Used for kruskal(), kruskalWithComponents(), prim() and primWithParents()
        void quickSortGraphs(Edge edge[], int left, int right)
        {
            if (left < right)
            {
                int middle = (left + right) / 2, i = left, j = right, d = 0;
                Edge aux = edge[left];
                edge[left] = edge[middle];
                edge[middle] = aux;
                while (i < j)
                {
                    if (edge[i].c > edge[j].c || (edge[i].c == edge[j].c && (edge[i].x > edge[j].x || (edge[i].x == edge[j].x && edge[i].y > edge[j].y))))
                    {
                        aux = edge[i];
                        edge[i] = edge[j];
                        edge[j] = aux;
                        d = 1 - d;
                    }
                    i += d;
                    j -= 1 - d;
                }
                quickSortGraphs(edge, left, i - 1);
                quickSortGraphs(edge, i + 1, right);
            }
        }

        // Used for royFloydWithPathShow()
        template <typename arr> inline
        void royFloydFindPath(int firstVertex, int secondVertex, int nextVertex[][MAX], int &length, arr &path)
        {
            length = 0;
            if (nextVertex[firstVertex][secondVertex])
            {
                path[0] = firstVertex;
                while (path[length] != secondVertex)
                {
                    path[length + 1] = nextVertex[path[length]][secondVertex];
                    ++length;
                }
                path[length++] = secondVertex;
            }
        }
    public :
        // O(nrVertices^2)
        template <typename arr> inline
        void dijkstra(int nrVertices, int startingVertex, int graph[][101], arr &distance, arr &selected)
        {
            //graph[i][i] = 0, graph[i][j] = INF, selected[i] = false before input
            //graph[x][y] = c during input
            int vertex, minimum, i;
            for (i = 1; i <= nrVertices; ++i)
                distance[i] = graph[startingVertex][i];
            selected[startingVertex] = true;
            do
            {
                vertex = -1;
                minimum = INFINITY;
                for (i = 1; i <= nrVertices; ++i)
                    if (!selected[i] && distance[i] < minimum)
                    {
                        minimum = distance[i];
                        vertex = i;
                    }
                if (vertex == -1)
                    break;
                selected[vertex] = true;
                for (i = 1; i <= nrVertices; ++i)
                    if (!selected[i] && distance[i] > distance[vertex] + graph[vertex][i])
                        distance[i] = distance[vertex] + graph[vertex][i];
            }
            while (1);
            for (i = 1; i <= nrVertices; ++i)
                if (distance[i] == INFINITY)
                    distance[i] = -1;
        }

        // O(nrVertices^2)
        template <typename arr> inline
        void dijkstraWithPath(int nrVertices, int startingVertex, int graph[][MAX], arr &distance, arr &previousVertex, arr &selected)
        {
            //graph[i][i] = 0, graph[i][j] = INFINITY, selected[i] = false before input
            //graph[x][y] = c during input
            int vertex, minimum, i;
            for (i = 1; i <= nrVertices; ++i)
            {
                if (graph[startingVertex][i] < INFINITY)
                    previousVertex[i] = startingVertex;
                distance[i] = graph[startingVertex][i];
            }
            selected[startingVertex] = true;
            do
            {
                vertex = -1;
                minimum = INFINITY;
                for (i = 1; i <= nrVertices; ++i)
                    if (!selected[i] && distance[i] < minimum)
                    {
                        minimum = distance[i];
                        vertex = i;
                    }
                if (vertex == -1)
                    break;
                selected[vertex] = true;
                for (i = 1; i <= nrVertices; ++i)
                    if (!selected[i] && distance[i] > distance[vertex] + graph[vertex][i])
                    {
                        previousVertex[i] = vertex;
                        distance[i] = distance[vertex] + graph[vertex][i];
                    }
            }
            while (1);
            for (i = 1; i <= nrVertices; ++i)
                if (distance[i] == INFINITY)
                    distance[i] = -1;
        }

        // Used with dijkstraWithPath()
        template <typename arr> inline
        void dijkstraWithPathShow(int nrVertices, int startingVertex, arr &distance, arr &path)
        {
            int i;
            for (i = 1; i <= nrVertices; ++i)
            {
                std::cout << distance[i];
                if (distance[i] != -1)
                {
                    int j, length;
                    dijkstraFindPath(i, startingVertex, distance, length, path);
                    std::cout << " : ";
                    for (j = 0; j < length; ++j)
                        std::cout << path[j] << ' ';
                }
                std::cout << '\n';
            }
        }

        // O(nrEdges * log(nrVertices))
        int kruskal(int nrVertices, int nrEdges, Edge edge[])
        {
            int i, j, nr = 0, cost = 0, aux, cc[MAX];
            for (i = 1; i <= nrVertices; ++i)
                cc[i] = i;
            quickSortGraphs(edge, 0, nrEdges - 1);
            for (i = 0; i < nrEdges && nr < nrVertices - 1; ++i)
                if (cc[edge[i].x] != cc[edge[i].y])
                {
                    cost += edge[i].c;
                    ++nr;
                    aux = cc[edge[i].y];
                    for (j = 1; j <= nrVertices; ++j)
                        if (cc[j] == aux)
                            cc[j] = cc[edge[i].x];
                }
            return cost;
        }

        // O(nrEdges * log(nrVertices))
        template <typename arr> inline
        void kruskalWithComponents(int nrVertices, int nrEdges, Edge edge[], Edge component[])
        {
            int i, j, nr = 0, cost = 0, aux, cc[MAX];
            for (i = 1; i <= nrVertices; ++i)
                cc[i] = i;
            quickSortGraphs(edge, 0, nrEdges - 1);
            for (i = 0; i < nrEdges && nr < nrVertices - 1; ++i)
                if (cc[edge[i].x] != cc[edge[i].y])
                {
                    cost += edge[i].c;
                    component[nr++] = edge[i];
                    aux = cc[edge[i].y];
                    for (j = 1; j <= nrVertices; ++j)
                        if (cc[j] == aux)
                            cc[j] = cc[edge[i].x];
                }
            std::cout << cost << '\n';
            for (i = 0; i < nrVertices - 1; ++i)
                std::cout << component[i].x << ' ' << component[i].y << '\n';
        }

        // O((nrVertices + nrEdges) * log(nrVertices))
        template <typename arr> inline
        int prim(int nrVertices, int nrEdges, int parentVertex, Edge edge[], arr &selected)
        {
            int cost = 0, i, j;
            for (i = 1; i <= nrVertices; ++i)
                selected[i] = false;
            selected[parentVertex] = true;
            quickSortGraphs(edge, 0, nrEdges - 1);
            for (i = 1; i < nrVertices; ++i)
            {
                j = 0;
                while (selected[edge[j].x] == selected[edge[j].y] && j < nrEdges)
                    ++j;
                cost += edge[j].c;
                if (!selected[edge[j].x])
                    selected[edge[j].x] = true;
                else
                    selected[edge[j].y] = true;
            }
            return cost;
        }

        // O((nrVertices + nrEdges) * log(nrVertices))
        template <typename arr> inline
        void primWithParentsShow(int nrVertices, int nrEdges, Edge edge[], arr &parent, arr &selected)
        {
            int i, j, cost;
            for (i = 1; i <= nrVertices; ++i)
            {
                primWithParents(nrVertices, nrEdges, i, cost, edge, parent, selected);
                if (i == 1)
                    std::cout << cost << '\n';
                std::cout << i << " : ";
                for (j = 1; j <= nrVertices; ++j)
                    std::cout << parent[j] << ' ';
                std::cout << '\n';
            }
        }

        // O(nrVertices^3)
        void royWarshall(int nrVertices, int adjacencyMatrix[][MAX])
        {
            //adjacencyMatrix[i][j] = 0 before input
            //adjacencyMatrix[x][y] = 1 during input
            int i, j, k;
            for (k = 1; k <= nrVertices; ++k)
                for (i = 1; i <= nrVertices; ++i)
                    for (j = 1; j <= nrVertices; ++j)
                        if (adjacencyMatrix[i][k] * adjacencyMatrix[k][j])
                            adjacencyMatrix[i][j] = 1;
        }

        // O(nrVertices^3)
        void royFloyd(int nrVertices, int graph[][MAX])
        {
            //graph[i][i] = 0, graph[i][j] = INFINITY before input
            //graph[x][y] = c during input
            int i, j, k;
            for (k = 1; k <= nrVertices; ++k)
                for (i = 1; i <= nrVertices; ++i)
                    for (j = 1; j <= nrVertices; ++j)
                        if (graph[i][j] > graph[i][k] + graph[k][j])
                            graph[i][j] = graph[i][k] + graph[k][j];
            for (int i = 1; i <= nrVertices; ++i)
                for (int j = 1; j <= nrVertices; ++j)
                    if (graph[i][j] == INFINITY)
                        graph[i][j] = -1;
        }

        // O(nrVertices^3)
        void royFloydWithPath(int nrVertices, int graph[][MAX], int nextVertex[][MAX])
        {
            //graph[i][i] = 0, graph[i][j] = INFINITY, nextVertex[i][j] = 0 before input
            //graph[x][y] = c, nextVertex[x][y] = y during input
            int i, j, k;
            for (k = 1; k <= nrVertices; ++k)
                for (i = 1; i <= nrVertices; ++i)
                    for (j = 1; j <= nrVertices; ++j)
                        if (graph[i][j] > graph[i][k] + graph[k][j])
                        {
                            graph[i][j] = graph[i][k] + graph[k][j];
                            nextVertex[i][j] = nextVertex[i][k];
                        }
        }

        // Used with royFloydWithPath
        template <typename arr> inline
        void royFloydWithPathShow(int nrVertices, arr &path, int nextVertex[][MAX])
        {
            int i, j, k, length;
            for (i = 1; i <= nrVertices; ++i)
                for (j = 1; j <= nrVertices; ++j)
                {
                    std::cout << i << ' ' << j;
                    royFloydFindPath(i, j, nextVertex, length, path);
                    if (length)
                    {
                        std::cout << " : ";
                        for (k = 0; k < length; ++k)
                            std::cout << path[k] << ' ';
                    }
                    std::cout << '\n';
                }
        }

        // O(nrVertices + nrEdges)
        template <typename arr> inline
        void topologicalSort(int nrVertices, int nrEdges, Edge edge[], arr &indegree, arr &queue)
        {
            // indegree[i] = 0 before input
            // ++indegree[edge[i].y] during input
            int i, vertex, left = 1, right = 1;
            for (i = 1; i <= nrVertices; ++i)
                if (!indegree[i])
                    queue[right++] = i;
            while (left <= right)
            {
                vertex = queue[left++];
                for (i = 1; i <= nrEdges; ++i)
                    if (edge[i].x == vertex)
                    {
                        --indegree[edge[i].y];
                        if (!indegree[edge[i].y])
                            queue[right++] = edge[i].y;
                    }
            }
        }
};


class numbersXXL_ALGO
{
    private :
        #define MAXLENGTH 1001
    public :
        template <typename arr> inline
        void readXXLNumber(int &numberLength, arr &numberInt)
        {
            int i;
            char numberChar[MAXLENGTH];
            std::cin >> numberChar;
            numberLength = strlen(numberChar);
            for (i = strlen(numberChar) - 1; i >= 0; --i)
                numberInt[numberLength - i - 1] = numberChar[i] - '0';
        }

        template <typename arr> inline
        void showXXLNumber(int numberLength, arr &number)
        {
            for (int i = numberLength - 1; i >= 0; --i)
                std::cout << number[i];
        }

        template <typename arr> inline
        void sumOfXXLNumbers(int firstLength, arr &firstNumber, int secondLength, arr &secondNumber, int &sumLength, arr &sum)
        {
            int i, transport = 0;
            if (firstLength < secondLength)
            {
                for (i = firstLength; i < secondLength; ++i)
                    firstNumber[i] = 0;
                sumLength = secondLength;
            }
            else
            {
                for (i = secondLength; i < firstLength; ++i)
                    secondNumber[i] = 0;
                sumLength = firstLength;
            }
            for (i = 0; i < sumLength; ++i)
                if (firstNumber[i] + secondNumber[i] + transport >= 10)
                {
                    sum[i] = firstNumber[i] + secondNumber[i] + transport - 10;
                    transport = 1;
                }
                else
                {
                    sum[i] = firstNumber[i] + secondNumber[i] + transport;
                    transport = 0;
                }
            if (transport)
                sum[sumLength++] = 1;
        }

        template <typename arr> inline
        void differenceOfXXLNumbers(int firstLength, arr &firstNumber, int secondLength, arr &secondNumber, int &differenceLength, arr &difference)
        {
            int i, transport = 0;
            for (i = secondLength; i < firstLength; ++i)
                secondNumber[i] = 0;
            differenceLength = firstLength;
            for (i = 0; i < differenceLength; ++i)
                if (firstNumber[i] - secondNumber[i] - transport < 0)
                {
                    difference[i] = firstNumber[i] - secondNumber[i] - transport + 10;
                    transport = 1;
                }
                else
                {
                    difference[i] = firstNumber[i] - secondNumber[i] - transport;
                    transport = 0;
                }
                while (!difference[differenceLength - 1] && differenceLength != 1)
                    --differenceLength;
        }

        template <typename arr> inline
        void productOfXXLNumbers(int firstLength, arr &firstNumber, int secondLength, arr &secondNumber, int &productLength, arr &product)
        {
            int i, j, transport = 0;
            productLength = secondLength + firstLength - 1;
            for (i = 0; i < productLength; ++i)
                product[i] = 0;
            for (i = 0; i < firstLength; ++i)
                for (j = 0; j < secondLength; ++j)
                    product[i + j] += firstNumber[i] * secondNumber[j];
            for (i = 0; i < productLength; ++i)
                product[i] += transport, transport = product[i] / 10, product[i] %= 10;
            if (transport)
                product[productLength++] = transport;
        }

        template <typename arr> inline
        void divisionOfXXLNumber(int &numberLength, arr &number, int div, int &dividendLength, arr &dividend)
        {
            dividendLength = numberLength;
            int i, transport = 0;
            for (i = numberLength - 1; i >= 0; --i)
            {
                dividend[i] = (transport * 10 + number[i]) / div;
                transport = (transport * 10 + number[i]) % div;
            }
            while (!dividend[dividendLength - 1])
                --dividendLength;
        }

        template <typename arr> inline
        void power2XXLNumber(int &power2Length, arr &power2, int power)
        {
            power2Length = 1, power2[0] = 1;
            int i, transport, x;
            while (power--)
            {
                transport = 0;
                for (i = 0; i < power2Length; ++i)
                {
                    x = transport + power2[i] * 2;
                    power2[i] = x % 10;
                    transport = x / 10;
                }
                if (transport)
                    power2[power2Length++] = 1;
            }
        }

        template <typename arr> inline
        void XXLNumberAsSumOfPow2(int numberLength, arr &number, int &exponentLength, arr &exponent)
        {
            exponentLength = 0;
            int nrOfBits = 0, i, transport, dv2[MAXLENGTH];
            do
            {
                transport = 0;
                if (number[0] % 2)
                    exponent[exponentLength++] = nrOfBits;
                ++nrOfBits;
                for (i = numberLength - 1; i >= 0; --i)
                {
                    dv2[i] = (transport * 10 + number[i]) / 2;
                    transport = (transport * 10 + number[i]) % 2;
                }
                if (!dv2[numberLength - 1])
                    --numberLength;
                for (i = 0; i < numberLength; ++i)
                    number[i] = dv2[i];
            } while (numberLength);
        }
};


class bitOperations_ALGO
{
    public :
        void setTrueBitAtPosition(int &number, int position)
        {
            number |= (1 << position);
        }

        bool checkIfTrueBitAtPosition(int number, int position)
        {
            return number & (1 << position);
        }

        void removeLastIndexTrueBits(int number, int index)
        {
            number &= (-1 << index);
        }

        void setFalseBitAtPosition(int &number, int position)
        {
            number &= ~(1LL << position);
        }
};


class benchmark_ALGO
{
    private :
        #define limit (int)1e5
        #define bm(x) benchmarch(x, #x)
    public :
        void benchmarch(void (*f)(), const char *name)
        {
            std::cout << "Running " << name << '\n';
            clock_t time = clock();
            f();
            time = clock() - time;
            std::cout << "Time taken = " << (double)time / CLOCKS_PER_SEC << "\n\n";
        }

        template <typename arr> inline
        void fillWithRandomNumbers(arr &array, int size)
        {
            srand(time(NULL));
            int left = limit, right = limit;
            for (int i = 0; i < size; ++i)
                array[i] = rand() % (left + right + 1) - rand() % left;
        }
};


#endif