---
Title: "Interaction graph"
weight: 3
---

The code can generate interaction graphs, which can help debugging, if results are not as expected. This can be enabled with the {{< variable "Debug" >}} variable set to {{< code "interaction_graph" >}}.

The following input file 

{{% expand "input file" %}}
{{< code-block >}}
#include_input testsuite/multisystem/02-interaction_graph.01-three_body.inp
{{< /code-block >}}
{{% /expand %}}

gives rise to the following graph:
{{% graphviz-file "/static/graph_data/interaction_graph.dot" %}}

This graph can be created using the {{< code "graphviz" >}} package.

{{% expand "How to install and use graphviz" %}}
On {{% name Debian %}}-like systems, {{< code graphviz >}} can simply be installed with
{{% code-block %}}
sudo apt-get install graphviz
{{% /code-block %}}
Then, the {{< code "interaction_graph.dot" >}} can be converted into a {{< code png >}} file by using:
{{% code-block %}}
dot -Tpng interaction_graph.dot > interaction_graph.png
{{% /code-block %}}

{{% /expand %}}