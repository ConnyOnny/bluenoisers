extern crate bluenoisers;
use bluenoisers::blue_noise;

fn main() {
    // create the samples
    let samples = blue_noise(vec![640.,480.], 4., 30);
    // create a boolean image to put the samples as pixels
    let mut image = vec![vec![false; 640]; 480];
    for s in samples {
        image[s[1] as usize][s[0] as usize] = true;
    }
    // write the image in portable bitmap format to stdout
    println!("P1");
    println!("{} {}", 640, 480);
    for y in 0..480 {
        for x in 0..640 {
            print!("{} ", if image[y][x] { "1" } else { "0" });
        }
        println!("");
    }
}
