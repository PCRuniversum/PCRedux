## PCRedux app

PCRedux app is a web server for extracting machine learning features from PCR curves using the [PCRedux](https://cran.r-project.org/package=PCRedux) package.

### Data format

The expected input data is <b>.csv</b> spreadsheet in the [qpcR](https://cran.r-project.org/package=qpcR) format. The first column represents the cycle number and the following commands represent particular runs

<table border="1" style="width:100%">
  <tr>
    <td><b>Cycle</b></td>
    <td><b>Run1</b></td> 
    <td><b>Run2</b></td> 
    <td>...</td> 
  </tr>
  <tr>
    <td>1</td>
    <td>17.22</td>
    <td>15.52</td>
    <td>...</td> 
  </tr>
  <tr>
    <td>2</td>
    <td>17.31</td>
    <td>17.17</td>
    <td>...</td> 
  </tr>
  <tr>
    <td>3</td>
    <td>16.99</td>
    <td>17.01</td>
    <td>...</td> 
  </tr>
  <tr>
    <td>...</td>
    <td>...</td>
    <td>...</td>
    <td>...</td> 
  </tr>
</table>
<br>

### Authors

- [Stefan Rödiger](https://www.researchgate.net/profile/Stefan_Roediger) | Stefan.Roediger[at]b-tu.de
- [Michał Burdukiewicz](https://www.researchgate.net/profile/Michal_Burdukiewicz) | michalburdukiewicz[at]gmail.com 
- [Dominik Rafacz](https://github.com/DominikRafacz) dominikrafacz[at]gmail.com
