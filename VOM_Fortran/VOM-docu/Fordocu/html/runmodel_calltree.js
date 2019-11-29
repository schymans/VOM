function runmodel(node, onCompleteCallback)
{
   var myobj = { label: "transpmodel", id: "transpmodel", href: "transpmodel.html", target:"basefrm" };
   var tmpNode21 = new YAHOO.widget.TextNode(myobj, node, false);
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}
