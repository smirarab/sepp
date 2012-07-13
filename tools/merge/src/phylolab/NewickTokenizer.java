package phylolab;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class NewickTokenizer
{
  String tree;
  Matcher pattern;
  boolean strip;

  public NewickTokenizer(String input)
  {
    init(input, true);
  }

  public NewickTokenizer(String input, boolean strip) {
    init(input, strip);
  }

  private void init(String input, boolean strip) {
    this.tree = input;
    this.strip = strip;
    if (strip)
      this.pattern = Pattern.compile("([(])|([)][^,:;)]*)|([;])|(:)|([^,);(:]*)").matcher(this.tree);
    else {
      this.pattern = Pattern.compile("([(])|([)][^,;)]*)|([;])|([^,);(]*)").matcher(this.tree);
    }
    this.pattern.find();
  }

  public boolean hasNext() {
    return !this.pattern.hitEnd();
  }
  public String nextToken() {
    String res = this.pattern.group();
    this.pattern.find();

    if ((this.strip) && (res.startsWith(")"))) {
      return ")";
    }
    if ("".equals(res))
      return nextToken();
    if (":".equals(res))
    {
      nextToken();
      return nextToken();
    }
    return res;
  }

  public static void main(String[] args)
  {
    String tree = "(2sd:32{2},((M:29{23},(A:10[2],(C:1,D:3)0.30[12]:232,B:12),N)22:232[12],L));";
    NewickTokenizer tokenizer = new NewickTokenizer(tree, false);
    while (tokenizer.hasNext())
      System.out.print(tokenizer.nextToken() + "  ");
  }
}