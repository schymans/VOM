function transpmodel_init(node, onCompleteCallback)
{
   var myobj = { label: "vom_init_vegpar", id: "vom_init_vegpar", href: "vom_init_vegpar.html", target:"basefrm" };
   var tmpNode23 = new YAHOO.widget.TextNode(myobj, node, false);
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}
